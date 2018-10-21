import ROOT
import math
import array


ROOT.gInterpreter.ProcessLine('''
const Int_t palette_size = 512;
Int_t MyPalette[palette_size];
Double_t Red[]    = {231/255.0,255/255.0,255/255.0,0/255.0,0/255.0,118/255.0};
Double_t Green[]  = {0/255.0,140/255.0,239/255.0,129/255.0,68/255.0,0/255.0};
Double_t Blue[]   = {0/255.0,0/255.0,0/255.0,31/255.0,255/255.0,137/255.0};
Double_t Length[] = {0,0.2,0.4,0.6,0.8,1.0};
Int_t FI = TColor::CreateGradientColorTable(6, Length, Red, Green, Blue, palette_size);
for(int i=0;i<palette_size;i++){
  MyPalette[i] = FI+(palette_size-1)-i;
}
gStyle->SetPalette(palette_size, MyPalette);
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
  std::vector<int> mData;
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
    FastHist self = FastHist(h->GetSize()-2, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetSize()-1));
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
  int GetBinContent(int i) const{
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
    std::function<float(int)> model = with_signal ? std::function<float(int)>([m, this](int i){return mModel.get_H1(m,i);}) : std::function<float(int)>([m, this](int i){return mModel.get_H0(i);});
    for(int i = 0; i < mModel.size();i++){
      int k = gRandom->Poisson(model(i));
      LL_H1 += log_likelihood_bin(mModel.get_H1(m, i), k);
      LL_H0 += log_likelihood_bin(mModel.get_H0(i), k);
      SetBinContent(i+1,k);
    }
    return std::vector({(float)mEntries, LL_H1, LL_H0});
  }
  std::vector<float> ComputeH1H0(float m = 2.1){
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
  std::vector<float> ComputeH1VariableMass(float m0, float m1, float m_delta, bool print = false){
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
    for(int& x: mData){
      x = 0;
    }
  }


  std::vector<float> do_psuedo_mt(int n_entries, bool with_Z, float m, bool dynamic = false) const{
    std::vector<float> ret(n_entries);
    int n_threads = std::thread::hardware_concurrency();
    std::thread thread[n_threads];
    for(int thread_id = 0; thread_id < n_threads; thread_id++){
      thread[thread_id] = std::thread(
        [&, thread_id]{
          FastHist my_h = *this;
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
};
''')

ROOT.gRandom.SetSeed(314158)
# ROOT.gInterpreter.ProcessLine('#include "FastHist.cpp"')
canv = ROOT.TCanvas("canv","plots for SDA course - ROOT 5", 2000, 1000)
fnc_bg = ROOT.TF1("background","[0]*exp(-x)/(exp(-1)-exp(-3))",1,3)
fnc_signal = ROOT.TF1("signal","[2]*TMath::Gaus(x, [0], [1],kTRUE)",1,3)
fnc_combined = ROOT.TF1("combined", "[2]*exp(-x)/(exp(-1)-exp(-3)) + [3]*TMath::Gaus(x, [0], [1],kTRUE)", 1,3)


n_bg = 200
n_signal = 10
h = ROOT.FastHist(60,1,3)

fnc_bg.SetParameter(0,n_bg*h.GetBinWidth(1));
#set m = 2.1
fnc_combined.SetParameter(0,2.1);
#set sigma
fnc_combined.SetParameter(1,50e-3);
#200 out of 210 are bg
fnc_combined.SetParameter(2,n_bg*h.GetBinWidth(1));
#10 out of 210 are signal
fnc_combined.SetParameter(3,n_signal*h.GetBinWidth(1));

def generate_psuedos():
    canv.Divide(2,1)
    canv.cd(1)
    h.FillRandom(True)
    root_h_H1 = ROOT.TH1F("Psuedo-experiment H1", "H1;invariant mass (TeV);entries per bin", 60, 1,3)
    h.copy_to_root_hist(root_h_H1)
    root_h_H1.Draw()
    fnc_combined.Draw("CSAME")
    fnc_bg.SetLineStyle(2)
    fnc_bg.Draw("CSAME")
    canv.cd(2)
    h.FillRandom(False)
    root_h_H0= ROOT.TH1F("Psuedo-experiment H0", "H0;invariant mass (TeV);entries per bin", 60, 1,3)
    h.copy_to_root_hist(root_h_H0)
    root_h_H0.Draw()
    fnc_bg.Draw("CSAME")
    canv.Update()
    canv.SaveAs("results/psuedo.png")
    canv.Clear()

def do_psuedo(with_Z, m = 2.1, dynamic = False):
    (n,H1,H0) = h.FillRandom(with_Z, m)
    if dynamic:
        (H1, m) = h.ComputeH1VariableMass(1.2, 2.8, 0.01)
    return H1-H0

def generate_single_LLR_distribution(with_particle, n_entries, m = 2.1, dynamic = False, hist_entries = None, hist_range = [-3,12]):
    if not hist_entries: hist_entries = n_entries/100
    LLR_H = ROOT.FastHist(hist_entries,*hist_range)
    results = h.do_psuedo_mt(n_entries, with_particle, m, dynamic)
    for x in results:
        LLR_H.Fill(x)
    return LLR_H

def generate_LLR_H1(n_entries, m = 2.1, dynamic = False, hist_entries = None, hist_range = [-3,12]):
    return generate_single_LLR_distribution(True, n_entries, m, dynamic, hist_entries, hist_range)

def generate_LLR_H0(n_entries, m = 2.1, dynamic = False, hist_entries = None, hist_range = [-3,12]):
    return generate_single_LLR_distribution(False, n_entries, m, dynamic, hist_entries, hist_range)

def generate_LLR_distribution(n_entries, m = 2.1, dynamic = False, hist_entries = None, hist_range = [-3,12]):
    LLR_H1 = generate_LLR_H1(n_entries, m, dynamic, hist_entries, hist_range)
    LLR_H0 = generate_LLR_H0(n_entries, m, dynamic, hist_entries, hist_range)
    return (LLR_H1, LLR_H0)

def plot_LLR_distribution(H1, H0, title):
    h_H0 = ROOT.TH1F("H0"+title, "Log likelihood-ratio;\lambda;Probablity", 1,0,1)
    h_H1 = ROOT.TH1F("H1"+title, "llr", 1,0,1)
    H0.copy_to_root_hist(h_H0)
    h_H0.Scale(1/(H0.GetEntries()*H0.GetBinWidth()))
    H1.copy_to_root_hist(h_H1)
    h_H1.Scale(1/(H1.GetEntries()*H1.GetBinWidth()))
    h_H0.SetLineColor(ROOT.kBlue)
    h_H1.SetLineColor(ROOT.kRed)
    h_H0.SetLineWidth(3)
    h_H1.SetLineWidth(3)
    h_H0.Draw("HIST")
    h_H1.Draw("HIST SAME")
    canv.Update()
    canv.SaveAs("results/"+title+".png")
    canv.Clear()

def read_real_data():
    h = ROOT.FastHist(60,1,3)
    for line in open("./assignment6-dataset.txt")  :   # loop over every line in the input file
        h.Fill(float(line))
    h.create_models(n_bg, n_signal)
    return h

#make a cdf of H0
def make_cdf(title, H0):
    cdf = ROOT.TH1F("cdf"+title, ";H0;Cumulative distribution", 1,0,1)
    H0.make_cdf().copy_to_root_hist(cdf)
    cdf.Scale(1.0/H0.GetEntries())
    cdf.SetLineWidth(3)
    cdf.Draw("HISTLF2")
    canv.Update()
    canv.SaveAs("results/"+title+".png")
    canv.Clear()
    return cdf

#read the real data
data = read_real_data()

#question a
generate_psuedos()
ROOT.gStyle.SetOptStat(0)           # Do not show statistics box

#question b
(H1, H0) = generate_LLR_distribution(1000,hist_range = [-10,14])
H1_mean = H1.GetMean()
print H1_mean, " = H1 mean" # 2.8868
print H0.GetMean(), " = H0 mean" # -2.3788
plot_LLR_distribution(H1, H0, "LLR")

#question c
cdf_H1 = make_cdf("cdf_H1", H1)
cdf = make_cdf("cdf", H0)
print 1.0 - cdf.Interpolate(H1_mean), " = predicted p-value" #0.0073
(H1, H0) = data.ComputeH1H0()
print 1.0 - cdf.Interpolate(H1-H0), " = measured p-value" # 0.00264800

#part two
(H1, H0) = generate_LLR_distribution(1000, 2.1, True)

H1_mean = H1.GetMean()
print H1_mean, " = H1 mean" # 2.8868
print H0.GetMean(), " = H0 mean" # -2.3788
plot_LLR_distribution(H1, H0, "LLR_dynamic")

cdf_dy = make_cdf("cdf_dynamic", H0)
cdf_H1_dy = make_cdf("cdf_H1_dy", H1)
cdf_dy.SetLineColor(ROOT.kRed)
cdf_H1_dy.SetLineColor(ROOT.kRed)
canv.Divide(2,1)
canv.cd(1)
cdf.Draw("HISTL")
cdf.SetTitle("CDF H0;H0;Cumulative distribution")
cdf_dy.Draw("HISTL SAME")
canv.cd(2)
cdf_H1.SetTitle("CDF H1;H1;Cumulative distribution")
cdf_H1.Draw("HISTL")
cdf_H1_dy.Draw("HISTL SAME")
canv.Update()
canv.SaveAs("results/cdfs_compared.png")
canv.Clear()
print 1.0 - cdf_dy.Interpolate(H1_mean), " = predicted dynamic p-value" #0.0073
H0 = data.ComputeH0()
(H1, m_fit) = data.ComputeH1VariableMass(1.2, 2.8, 0.1)
print 1.0 - cdf_dy.Interpolate(H1-H0), " = measured dynamic p-value" # 0.00264800
print m_fit, " = m fit"

fnc_combined.SetParameter(0,m_fit)
root_h = ROOT.TH1F("fitted m", "Fit result;invariant mass (TeV);entries per bin", 60, 1,3)
data.copy_to_root_hist(root_h)
root_h.Draw()
fnc_combined.Draw("CSAME")
fnc_bg.SetLineStyle(2)
fnc_bg.Draw("CSAME")
canv.Update()
canv.SaveAs("results/fit.png")
canv.Clear()

nM = 10
M = array.array('d',[1.2 + 1.6*m_i/nM for m_i in range(nM + 1)])
P = array.array('d',[0 for _ in M])
n_entries = 500
H1_2d = ROOT.TH2F('2dH1', 'LLR distribution of H1 vs particle mass M_{z};\lambda;M_{z};number of entries', 2000, -2, 16, nM, 1.2, 2.8)
for i,m in enumerate(M):
    H1 = generate_LLR_H1(n_entries, m, True, 2000, [-2,16])
    for j in range(n_entries):
        H1_2d.SetBinContent(j+1, i+1, H1.GetBinContent(j+1))
    H1_mean = H1.GetMean()
    p =  1.0 - cdf_dy.Interpolate(H1_mean)
    P[i] = p
    print m, H1_mean, p

H1_2d.RebinX(20)
H1_2d.Draw("COLZ")
canv.Update()
canv.SaveAs("results/2dH1H0.png")
canv.Clear()

graph = ROOT.TGraphErrors(len(M),M,P)
graph.SetMarkerColor(ROOT.kRed)
graph.SetMarkerStyle(21)
graph.SetTitle("Expected p-value as a function of true M;True Invariant Mass (TeV);Expected p-value")
graph.Draw("ALP")
canv.Update()
canv.SaveAs("results/expected_p_vs_M.png")
canv.Clear()
