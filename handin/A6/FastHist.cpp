#include <vector>
#include <algorithm>
#include <thread>

//I'm sorry, had to try. Became a little bit of a mess, but it works :)
//A class for caching the model predictions  so that the expensive TMath calls
// are not repeatedly evaluated.
class Model{
  std::vector<float> mH0Cache;
  std::map<float,std::vector<float>> mH1Cache;
public:
  Model(int n_bins) : mH0Cache(n_bins){}
  int size() const{ return mH0Cache.size();}
  void fill_bg(int n_bg, float x0, float binwidth){
    for(int x_i = 0; x_i < mH0Cache.size(); x_i++){
      float x = x0 + x_i*binwidth;
      mH0Cache[x_i] = binwidth*n_bg*TMath::Exp(-x)/(TMath::Exp(-1)-TMath::Exp(-3));
    }
  }
  void fill_model(int n_signal, float x0, float binwidth, float m){
    mH1Cache[m] = std::vector<float>(mH0Cache.size());
    for(int x_i = 0; x_i < mH0Cache.size(); x_i++){
      float x = x0 + (x_i+0.5)*binwidth;
      mH1Cache[m][x_i] = mH0Cache[x_i] + binwidth*n_signal*TMath::Gaus(x, m, 50e-3,kTRUE);
    }
  }
  float get_H1(float m, int i) const{
    return mH1Cache.at(m)[i];
  }
  float get_H0(int i) const{
    return mH0Cache[i];
  }
  bool contains_mass(float m) const{
    bool contains = mH1Cache.count(m) > 0;
    // if(!contains)
      // std::cout << "model contains " << m << "? " << contains << std::endl;
    return contains;
  }
};

//A simplified histogram class that exposes similar calls to ROOT
class FastHist{
  int mEntries = 0;
  //the entries
  std::vector<double> mData;
  //cache for model evaluation
  Model mModel;
  //bad hardcoding.
  int n_signal = 10;
  int n_bg = 200;

  float mStart = 0;
  float mDelta = 0;
  void fill_model_cache(int n_bg, int n_signal, float m = 2.1){
    mModel.fill_bg(n_bg, mStart, mDelta);
    mModel.fill_model(n_signal, mStart, mDelta, m);
  }
  void fill_bg_cache(int n_bg){
    mModel.fill_bg(n_bg, mStart, mDelta);
  }
  float log_likelihood_bin(float x, int k) const{
    return -x + k*log(x) - lgamma(k+1);
  }
public:
  FastHist(int n_bins = 0, float start = 0, float end = 1) : mData(n_bins+2, 0), mStart(start), mDelta((end-start)/n_bins), mModel(n_bins){
    mModel.fill_bg(n_bg, mStart, mDelta);
  }
  static FastHist FromROOT(const TH1F* h){
    auto self = FastHist(h->GetSize()-2, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetSize()-1));
    for(int i = 0; i < h->GetSize();i++){
      self.SetBinContent(i, h->GetBinContent(i));
    }
    return self;
  }
  void Fill(float x){
    mEntries += 1;
    int bin_id = (x - mStart + mDelta)/mDelta;
    bin_id = std::clamp(bin_id, 0, (int)(mData.size()-1));
    mData[bin_id] += 1;
  }
  auto GetBinContent(int i) const{
    return mData[i];
  }
  void SetBinContent(int i, int x){
    mEntries += x - mData[i];
    mData[i] = x;
  }
  int GetEntries() const{
    return mEntries;
  }
  int GetSize() const{
    return mData.size();
  }
  float GetBinWidth(int i = 0) const{
    return mDelta;
  }
  float GetBinCenter(int i) const{
    return mStart + mDelta*(i-0.5);
  }
  double GetMean() const{
    int total_entries = 0;
    double mean = 0;
    for(int i = 1; i < mData.size()-1; i++){
      total_entries += mData[i];
      mean += mData[i]*GetBinCenter(i);
    }
    return mean/total_entries;
  }
  //copies its contents to `h`.
  void copy_to_root_hist(TH1F *h){
    h->SetBins(GetSize()-2, mStart, mStart+(GetSize()-2)*mDelta);
    for(int i = 0; i < GetSize();i++)
        h->SetBinContent(i,GetBinContent(i));
    h->SetEntries(mEntries);
  }
  void create_models(int n_bg, int n_signal, bool with_signal = true, float m = 2.1){
    fill_model_cache(n_bg,n_signal, m);
  }

  //fills according to the model. With signal determines if a Z particle signal should be generated with mass m
  // returns the number of entries and the log likelihood of H1 and H0.
  std::vector<float> FillRandom(bool with_signal = true, float m = 2.1){
    float LL_H0 = 0.0;
    float LL_H1 = 0.0;
    if(!mModel.contains_mass(m)){
      mModel.fill_model(n_signal, mStart, mDelta, m);
    }
    auto model = with_signal ? std::function<float(int)>([m, this](int i){return mModel.get_H1(m,i);}) : std::function<float(int)>([m, this](int i){return mModel.get_H0(i);});
    for(int i = 0; i < mModel.size();i++){
      int k = gRandom->Poisson(model(i));
      LL_H1 += log_likelihood_bin(mModel.get_H1(m, i), k);
      LL_H0 += log_likelihood_bin(mModel.get_H0(i), k);
      SetBinContent(i+1,k);
    }
    return std::vector({(float)mEntries, LL_H1, LL_H0});
  }
  std::vector<float> ComputeH1H0(float m = 2.1) const{
    if(!mModel.contains_mass(m)){
      mModel.fill_model(n_signal, mStart, mDelta, m);
    }
    float LL_H0 = 0.0;
    float LL_H1 = 0.0;
    for(int i = 0; i < mModel.size();i++){
      LL_H1 += log_likelihood_bin(mModel.get_H1(m, i), mData[i+1]);
      LL_H0 += log_likelihood_bin(mModel.get_H0(i), mData[i+1]);
    }
    return std::vector({LL_H1,LL_H0});
  }

  float ComputeH1(float m = 2.1){
    if(!mModel.contains_mass(m)){
      mModel.fill_model(n_signal, mStart, mDelta, m);
    }
    float LL_H1 = 0.0;
    for(int i = 0; i < mModel.size();i++){
      LL_H1 += log_likelihood_bin(mModel.get_H1(m, i), mData[i+1]);
    }
    return LL_H1;
  }

  //computes H1 for m0:m1:m_delta and returns the best scoring log-likelihood corresponding m
  auto ComputeH1VariableMass(float m0, float m1, float m_delta, bool print = false){
    float best = -1e30;
    float best_m = 0;
    for(float m = m0; m <= m1; m+= m_delta){
      float LL_H1 = ComputeH1(m);
      if( LL_H1 > best ){
        best = LL_H1;
        best_m = m;
      }
    }
    return std::vector<float>({best, best_m});
  }
  float ComputeH0() const{
    float LL_H0 = 0.0;
    for(int i = 0; i < mModel.size();i++){
      LL_H0 += log_likelihood_bin(mModel.get_H0(i), mData[i+1]);
    }
    return LL_H0;
  }

  //Makes a non-normalized cdf of the current data.
  FastHist make_cdf() const{
    FastHist ret(mData.size()-2, mStart, mStart+mDelta*(mData.size()-2));
    float sum = 0;
    for(int i = 1; i <  mData.size()-1;i++){
      sum += mData[i];
      ret.SetBinContent(i,sum);
    }
    return ret;
  }
  //set everything to 0
  void Clear(){
    mEntries = 0;
    for(auto& x: mData){
      x = 0;
    }
  }

};

//Fills h with a psuedoexperiment according to the parameters given, and runs multi-threaded.
// returns a vector of lambdas.
std::vector<float> do_psuedo_mt(const FastHist *h, int n_entries, bool with_Z, float m, bool dynamic = false){
  std::vector<float> ret(n_entries);
  int n_threads = std::thread::hardware_concurrency();
  std::thread thread[n_threads];
  for(int thread_id = 0; thread_id < n_threads; thread_id++){
    thread[thread_id] = std::thread(
      [&, thread_id]{
        FastHist my_h = *h;
        for(int i = thread_id; i < n_entries; i += n_threads){
            auto fill_result = my_h.FillRandom(with_Z, m);
            float H1 = fill_result[1];
            float H0 = fill_result[2];
            if(dynamic)
              H1 = my_h.ComputeH1VariableMass(1.2, 2.8, 0.01)[0];
            ret[i] = H1-H0;
        }
      }
    );

  }
  for(int thread_id = 0; thread_id < n_threads; thread_id++){
    thread[thread_id].join();
  }
  return ret;
}
