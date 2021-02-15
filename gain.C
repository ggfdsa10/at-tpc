#include<stdio.h>
#include <iostream>
#include <fstream>


using namespace std;
using namespace TMath;

void gain()
{
  ofstream finaldata;
  finaldata.open("./output/Data_information.txt");

  int epoch = 67;
  
  double total_xmean = 0;
  double total_xsigma = 0;
  double total_ymean = 0;
  double total_ysigma = 0;
  double x_distence = 0;
  double y_distence = 0;

  TCanvas *c0 = new TCanvas("c0","Distence Distribution",1000,500);
  c0 -> Divide(2,1);

  TCanvas *cc = new TCanvas("cc","Electron Generate",600,500);
  cc -> Divide(1,1);

  TH1D *distence1 = new TH1D("X-Data","X-Distence Distribution",500, -3, 3);
  TH1D *distence2 = new TH1D("Y-Data","Y-Distence Distribution",500, -3, 3);
  TF1 *df1 = new TF1("df1","gaus",-1,1);
  TF1 *df2 = new TF1("df2","gaus",-1,1);

  TH1D *e1 = new TH1D("Data","Electron Generate Distribution", 100, 0,300000);
  

  for (int i =0; i<epoch; i++){
    double x0,y0,z0,x2,y2,z2,t2,e2;
    
    char name[epoch];
    sprintf(name, "./tpc_data/data_%d.txt", i+1);

    fstream file;
    file.open(name, ios::in);
    
    TCanvas *c1 = new TCanvas("c1", "distribution",1000,500);
    c1 -> Divide(2,1);
    
    TH1D *h1 = new TH1D("X-Data","X-distribution", 1300, -1,1);
    TH1D *h2 = new TH1D("Y-Data","Y-distribution", 1300, -1,1);

    TF1 *f1 = new TF1("f1", "gaus",-1,1);
    TF1 *f2 = new TF1("f2", "gaus",-1,1);

    int electron = 0;
    
    bool empty = file.peek() == EOF;
    if (empty == false){
      while (!file.eof()){
	file >> x0 >> y0 >> z0 >> x2 >> y2 >> z2 >> t2 >> e2;
	h1 -> Fill(x2);
	h2 -> Fill(y2);
	electron +=1;
      }
    }

    e1 -> Fill(electron);
    
    double xmean = h1 -> GetMean();
    double xsigma = h1 -> GetStdDev();
    double ymean = h2 -> GetMean();
    double ysigma = h2 -> GetStdDev();

    finaldata << " Epoch : " << i+1 << endl;
    finaldata << " X_Mean : " << xmean << endl;
    finaldata << " Y_Mean : " << ymean << endl;
    finaldata << " X_Standard Diviation : " << xsigma << endl;
    finaldata << " Y_Standard Diviation : " << ysigma << endl << endl;

    total_xmean += xmean;
    total_xsigma += xsigma;
    total_ymean += ymean;
    total_ysigma += ysigma;
    
    distence1 -> Fill(xmean-x0);
    distence2 -> Fill(ymean-y0);
    
    f1 -> SetParameter(0, 1/(sqrt(2*Pi())*xsigma));
    f1 -> SetParameter(1, xmean);
    f1 -> SetParameter(2, xsigma);

    f2 -> SetParameter(0, 1/(sqrt(2*Pi())*ysigma));
    f2 -> SetParameter(1, ymean);
    f2 -> SetParameter(2, ysigma);

    c1 -> cd(1);
    h1 -> GetXaxis() -> SetTitle("X [cm]");
    h1 -> GetYaxis() -> SetTitle("Number of electrons");
    h1 -> GetXaxis() -> SetRangeUser(xmean - 5*xsigma, xmean + 5*xsigma); 
    h1 -> SetFillColor(kBlue);
    h1 -> Draw();
    h1 -> Fit("f1","R");
    
    c1 -> cd(2);
    h2 -> GetXaxis() -> SetTitle("Y [cm]");
    h2 -> GetYaxis() -> SetTitle("Number of electrons");
    h2 -> GetXaxis() -> SetRangeUser(ymean - 5*ysigma, ymean + 5*ysigma);
    h2 -> SetFillColor(kBlue);
    h2 -> Draw();
    h2 -> Fit("f2","R");

    TLine *l1 = new TLine(x0, 0, x0, h1->GetMaximum());
    TLine *l2 = new TLine(y0, 0, y0, h2->GetMaximum());
    c1 -> cd(1);
    l1 -> SetLineWidth(2);
    l1 -> SetLineColor(kOrange);
    l1 -> Draw();
    
    c1 -> cd(2);
    l2 -> SetLineWidth(2);
    l2 -> SetLineColor(kOrange);
    l2 -> Draw();

    TLegend *leg1 = new TLegend(0.78,0.67,0.98,0.775);
    c1 ->cd(1);
    leg1 -> AddEntry(h1, "Electron Gain Data", "f");
    leg1 -> AddEntry(f1, "Fit Gaussian","l");
    leg1 -> AddEntry(l1, "Initial Position", "l");
    leg1 -> Draw();

    TLegend *leg2 = new TLegend(0.78,0.67,0.98,0.775);
    c1 ->cd(2);
    leg2 -> AddEntry(h1, "Electron Gain Data", "f");
    leg2 -> AddEntry(f1, "Fit Gaussian","l");
    leg2 -> AddEntry(l1, "Initial Position", "l");
    leg2 -> Draw();

    char xdistence[1];
    sprintf(xdistence, "Mean-Initial Distence = %f",abs(xmean-x0));
    TLatex *latex1 = new TLatex();
    c1 -> cd(1);
    if (h1 -> GetEntries() ==0){
      latex1 -> SetTextSize(0.1);
      latex1 -> DrawLatex(-0.6, 0.5, "Empty Data");
    }
    else{
      latex1 -> SetTextSize(0.02);
      latex1 -> DrawLatex(xmean,h1->GetMaximum(),xdistence);
    }


    char ydistence[1];
    sprintf(ydistence, "Mean-Initial Distence = %f",abs(ymean-y0));
    TLatex *latex2 = new TLatex();
    c1 -> cd(2);
    if (h1 -> GetEntries() ==0){
      latex2 -> SetTextSize(0.1);
      latex2 -> DrawLatex(-0.6, 0.5, "Empty Data");
    }
    else{
      latex2 -> SetTextSize(0.02);
      latex2 -> DrawLatex(ymean,h1->GetMaximum(),ydistence);
    }
    c1 -> Draw();

    char pic[epoch];
    sprintf(pic, "./output/data_%d.png", i+1);
    
    c1 ->SaveAs(pic,"PNG");
  }

  double x_mean_mean = distence1 -> GetMean();
  double x_mean_sigma = distence1 -> GetStdDev();
  double y_mean_mean = distence2 -> GetMean();
  double y_mean_sigma = distence2 -> GetStdDev();

  df1 -> SetParameter(0, 1/(sqrt(2*Pi())*x_mean_sigma));
  df1 -> SetParameter(1, x_mean_mean);
  df1 -> SetParameter(2, x_mean_sigma);

  df2 -> SetParameter(0, 1/(sqrt(2*Pi())*y_mean_sigma));
  df2 -> SetParameter(1, y_mean_mean);
  df2 -> SetParameter(2, y_mean_sigma);
  
  c0 -> cd(1);
  distence1 -> SetFillColor(kBlue);
  distence1 -> Draw();
  distence1 -> GetXaxis() -> SetTitle("X [cm]");
  distence1 -> GetYaxis() -> SetTitle("#lbar Mean - X_Initial #lbar Count");
  distence1 -> GetXaxis() -> SetRangeUser(x_mean_mean - 4*x_mean_sigma, x_mean_mean + 4*x_mean_sigma);
  distence1 -> Fit("df1","R");

  c0 ->cd(2);
  distence2 -> SetFillColor(kBlue);
  distence2 -> Draw();
  distence2 -> GetXaxis() -> SetTitle("Y [cm]");
  distence2 -> GetYaxis() -> SetTitle("#lbar Mean - Y_Initial #lbar Count");
  distence2 -> GetXaxis() -> SetRangeUser(y_mean_mean - 4*y_mean_sigma, y_mean_mean + 4*y_mean_sigma);
  distence2 -> Fit("df2","R");

  TLegend *dleg1 = new TLegend(0.78,0.67,0.98,0.775);
  c0 ->cd(1);
  dleg1 -> AddEntry(distence1, "Mean-X Distence", "f");
  dleg1 -> AddEntry(df1, "Fit Gaussian","l");
  dleg1 -> Draw();

  TLegend *dleg2 = new TLegend(0.78,0.67,0.98,0.775);
  c0 ->cd(2);
  dleg2 -> AddEntry(distence2, "Mean-Initial Distence", "f");
  dleg2 -> AddEntry(df2, "Fit Gaussian","l");
  dleg2 -> Draw();
  
  c0 -> Draw();
  c0 -> SaveAs("./output/Distence_mean.png","PNG");

  cc ->cd(1);
  e1 -> SetFillColor(kRed);
  e1 -> GetYaxis() -> SetRange(1,1000);
  e1 -> GetXaxis() -> SetTitle("Number of Electrons");
  e1 -> GetYaxis() -> SetTitle("Epoch Count");
  e1 -> Draw();
  cc -> Draw();
  gPad -> SetLogy();
  cc -> SaveAs("./output/Electron_Generate.png","PNG");
  
  cout << " Total Epoch : " << epoch << endl;
  cout << " Average X-Mean : " << total_xmean/epoch << endl;
  cout << " Average Y-Mean : " << total_ymean/epoch << endl;
  cout << " Average X-Standard Deviation : " << total_xsigma/epoch << endl;
  cout << " Average Y-Standard Deviation : " << total_ysigma/epoch << endl;

  finaldata << " Total Epoch : " << epoch << endl;
  finaldata << " Average X-Mean : " << total_xmean/epoch << endl;
  finaldata << " Average Y-Mean : " << total_ymean/epoch << endl;
  finaldata << " Average X-Standard Deviation : " << total_xsigma/epoch << endl;
  finaldata << " Average Y-Standard Deviation : " << total_ysigma/epoch << endl;

  finaldata.close();
}

  

  
