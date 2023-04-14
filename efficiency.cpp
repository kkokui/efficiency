#include <stdio.h>
#include <EdbDataSet.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TH1.h>

int main(int argc , char *argv[]){
	if(argc<=0){
		printf("Usage : ./efficiency\n");
		return 1;
	}

	EdbDataProc *dproc = new EdbDataProc;
	EdbPVRec *pvr = new EdbPVRec;
	// dproc->ReadTracksTree(*pvr, "linked_tracks.root",Form("nseg>=5&& abs(t.eTX + 0.01) < 0.01 && abs(t.eTY-0.004) < 0.01"));
	dproc->ReadTracksTree(*pvr, "linked_tracks.root", "nseg>=5");
	TEfficiency *pEff_angle =0;
	TEfficiency *pEff_plate =0;
	
	int ntrk = pvr->Ntracks();
	double x1, y1, z1, x2, y2, z2, angle;
	double bins[] = {0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.25,0.3,0.35,0.4,0.45,0.5};
	
	
	// TH1D *h_angle_total = new TH1D("total",";angle;",10,0,0.5);
	// TH1D *h_angle_passed = new TH1D("passed",";angle;",10,0,0.5);
	TH1D *h_angle_total = new TH1D("hist_angle_total", ";tan#theta;", 26, bins);
	TH1D *h_angle_passed = new TH1D("hist_angle_passed", ";tan#theta;", 26, bins);
	TH1D *h_plate_total = new TH1D("hist_plate_total", ";plate;", 90, 50,140);
	TH1D *h_plate_passed = new TH1D("hist_plate_passed", ";plate;", 90, 50,140);

	for(int itrk=0; itrk<ntrk; itrk++){
		EdbTrackP *t = pvr->GetTrack(itrk);
		int nseg = t->N();
		for(int iplate=50;iplate<=140;iplate++){
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
				if(s->Plate()==iplate) hitsOnThePlate=1;
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
			}
			
		}
		
	}
	TCanvas *c = new TCanvas();
	c->Print("hist_efficiency.pdf[");
	c->SetLogy(1);
	h_angle_total->Draw();
	c->Print("hist_efficiency.pdf");
	c->SetLogy(0);
	pEff_angle = new TEfficiency(*h_angle_passed, *h_angle_total);
	pEff_angle->SetTitle("Efficiency for each angle;tan#theta;efficiency");
	pEff_angle->Draw();
	c->Print("hist_efficiency.pdf");
	h_plate_total->Draw();
	c->Print("hist_efficiency.pdf");
	pEff_plate = new TEfficiency(*h_plate_passed, *h_plate_total);
	pEff_plate->SetTitle("Efficiency for each plate;plate;efficiency");
	pEff_plate->Draw();
	c->Print("hist_efficiency.pdf");
	c->Print("hist_efficiency.pdf]");
	TFile fout("efficiency.root","recreate");
	pEff_angle->Write();
	fout.Close();
	return 0;
}