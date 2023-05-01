#include <stdio.h>
#include <fstream>
#include <sstream>
#include <EdbDataSet.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TH1.h>
std::string module;
int zone, plMin, plMax, reco;
double Xcenter, Ycenter;
int read_volume_info()
{
	std::ifstream file;
	file.open("volume_info.txt");	
	if(file.fail())
	{
		printf("Failed to open volume_info.txt\n");
		return 0;
	}
	std::string line;
	int count=0;
	while (std::getline(file, line))
	{
		std::vector<std::string> elems;
		std::stringstream ss(line);
		std::string item;
		while(std::getline(ss, item, '='))
		{
			if(!item.empty())
				elems.push_back(item);
		}
		if("module"==elems.at(0))
		{
			module = elems.at(1);
			count++;
		}
		if ("zone" == elems.at(0))
		{
			zone = std::stoi(elems.at(1));
			count++;
		}
		if ("plMin" == elems.at(0))
		{
			plMin = std::stoi(elems.at(1));
			count++;
		}
		if ("plMax" == elems.at(0))
		{
			plMax = std::stoi(elems.at(1));
			count++;
		}
		if ("reco" == elems.at(0))
		{
			reco = std::stoi(elems.at(1));
			count++;
		}
		if ("Xcenter" == elems.at(0))
		{
			Xcenter = std::stod(elems.at(1));
			count++;
		}
		if ("Ycenter" == elems.at(0))
		{
			Ycenter = std::stod(elems.at(1));
			count++;
		}
		
	}
	// printf("%s,%d,%d,%d,%d,%f,%f\n",module.c_str(),zone,plMin,plMax,reco,Xcenter,Ycenter);
	if (7 == count)
	{
		printf("volume_info.txt read successfully\n");
		return 1;
	}
	else
	{
		return 0;
	}
}

int main(int argc , char *argv[]){
	if(argc<=2){
		printf("Usage : ./efficiency linked_tracks.root title \"nseg>=5\"\n");
		return 0;
	}
	if(0==read_volume_info())
	{
		printf("Could not read volume_info.txt\n");
		return 0;
	}
	TString linked_tracks = argv[1];
	TString title = argv[2];
	TString cut = argv[3];
	EdbDataProc *dproc = new EdbDataProc;
	EdbPVRec *pvr = new EdbPVRec;
	// dproc->ReadTracksTree(*pvr, "linked_tracks.root",Form("nseg>=5&& abs(t.eTX + 0.01) < 0.01 && abs(t.eTY-0.004) < 0.01"));
 	dproc->ReadTracksTree(*pvr, linked_tracks, cut);
	TEfficiency *pEff_angle =0;
	TEfficiency *pEff_plate =0;
	
	int ntrk = pvr->Ntracks();
	double x1, y1, z1, x2, y2, z2, angle;
	double bins[] = {0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.25,0.3,0.35,0.4,0.45,0.5};
	int nbins = 26;
	
	// TH1D *h_angle_total = new TH1D("total",";angle;",10,0,0.5);
	// TH1D *h_angle_passed = new TH1D("passed",";angle;",10,0,0.5);
	TH1D *h_angle_total = new TH1D("hist_angle_total", title+";tan#theta;", nbins, bins);
	TH1D *h_angle_passed = new TH1D("hist_angle_passed", title+";tan#theta;", nbins, bins);
	TH1D *h_plate_total = new TH1D("hist_plate_total", title+";plate;", plMax-plMin, plMin,plMax);
	TH1D *h_plate_passed = new TH1D("hist_plate_passed", title + ";plate;", plMax - plMin, plMin, plMax);
	TFile ntfout(Form("nt_not_passed_%s.root",title.Data()), "recreate");
	double x, y;
	TNtuple *nt = new TNtuple("nt","title","trackID:x:y:angle:pl:nseg");

	for(int itrk=0; itrk<ntrk; itrk++){
		EdbTrackP *t = pvr->GetTrack(itrk);
		int nseg = t->N();
		for(int iplate=plMin;iplate<=plMax;iplate++){
			int counts = 0;
			int hitsOnThePlate = 0;
			for(int iseg=0;iseg<nseg;iseg++){
				EdbSegP *s = t->GetSegment(iseg);
				
				if(s->Plate()==iplate-2||s->Plate()==iplate+2) counts++;
				if(s->Plate()==iplate-1){
					x1=s->X();
					y1=s->Y();
					z1=s->Z();
					counts++;
				}
				if(s->Plate()==iplate+1){
					x2=s->X();
					y2=s->Y();
					z2=s->Z();
					counts++;
				}
				if(s->Plate()==iplate)
				{
					hitsOnThePlate=1;
				}
			}
			if(counts==4) {
				angle = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))/(z2-z1);
				h_angle_total->Fill(angle);
				h_plate_total->Fill(iplate);
				if(hitsOnThePlate==1) 
				{
					h_angle_passed->Fill(angle);
					h_plate_passed->Fill(iplate);
				}
				else
				{
					nt->Fill(t->ID(),(x1+x2)/2,(y1+y2)/2,angle,iplate,nseg);
				}
			}
			
		}
		
	}
	nt->Write();
	TCanvas *c = new TCanvas();
	c->Print(Form("hist_efficiency_%s.pdf[", title.Data()));
	c->SetLogy(1);
	h_angle_total->Draw();
	c->Print(Form("hist_efficiency_%s.pdf",title.Data()));
	c->SetLogy(0);
	pEff_angle = new TEfficiency(*h_angle_passed, *h_angle_total);
	pEff_angle->SetTitle(Form("Efficiency for each angle (%s);tan#theta;efficiency",title.Data()));
	pEff_angle->Draw();
	c->Print(Form("hist_efficiency_%s.pdf",title.Data()));
	h_plate_total->Draw();
	c->Print(Form("hist_efficiency_%s.pdf",title.Data()));
	pEff_plate = new TEfficiency(*h_plate_passed, *h_plate_total);
	pEff_plate->SetTitle(Form("Efficiency for each plate (%s);tan#theta;efficiency", title.Data()));
	pEff_plate->Draw();
	c->Print(Form("hist_efficiency_%s.pdf",title.Data()));
	c->Print(Form("hist_efficiency_%s.pdf]", title.Data()));
	TFile fout(Form("efficiency_%s.root",title.Data()),"recreate");
	pEff_angle->Write();
	pEff_plate->Write();
	fout.Close();

	FILE *ftxt = fopen("efficiency.txt","w");
	for(int i=1;i<=nbins;i++)
	{
		fprintf(ftxt, "%.3f, ",(bins[i]+bins[i-1])/2);
	}
	fprintf(ftxt,"\n");
	for(int i=1;i<=nbins;i++)
	{
		fprintf(ftxt,"%f ",pEff_angle->GetEfficiency(i));
	}
	fprintf(ftxt,"\n");
	fclose(ftxt);
	return 1;
}