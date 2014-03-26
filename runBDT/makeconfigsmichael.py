#!/bin/python

import os, sys
import ConfigParser

def makeconfigfile(configfilename,NN_book_options,NN_vars,train_test,backgroundfilelist,signalsample):
    config = ConfigParser.RawConfigParser()
    
    config.add_section('NN')
    config.set('NN', "NN_book_options", NN_book_options)
    config.set('NN', "NN_vars", NN_vars)
    config.set('NN', "backgroundfilelist", backgroundfilelist)
    config.set('NN', "signalfilelist", signalsample)
    config.set('NN', "NN_train_test", train_test)
    
    with open(configfilename, 'wb') as configfile:
        config.write(configfile)

def getconfigname(counter,basedir):
    if not os.path.isdir(basedir):
        os.makedirs(basedir)
    return basedir + "setup_" + str(counter) + ".cfg"

if __name__ == "__main__":

    basedir = ""
    if not len(sys.argv) is 2:
        print "will use default directory for setups"
        basedir = os.getcwd() + "/setup/"
    else:
        basedir = sys.argv[1] #dir where you want the setups

    print basedir
   

    NN_vars_list = [
		    #t2bw	
                    #"met,lepton_pT,mlb_hemi,m3b,mT2W,b1_pt,dPhi_JetMet,dR_LepB" # laras stuff 100 
                    #"met,lepton_pT,mlb_hemi,m3b,mT2W,HT,dPhi_JetMet,jet1_pT" # laras stuff 101
                    #"met,lepton_pT,njets,mlb_hemi,m3b,mT2W,b1_pt,dPhi_JetMet,dR_LepB" # laras stuff 102 
                    #"met,lepton_pT,njets,mlb_hemi,m3b,mT2W,HT,dPhi_JetMet,jet1_pT" # laras stuff 103

		    #t2tt
		    #"met,mT2W,HTfrac,njets,dPhi_JetMet,lepton_pT,dR_LepB,Chi2SNT,jet1_pT"  #104
		    "met,lepton_pT,mlb_hemi,HTfrac,njets,dPhi_JetMet,b1_pt" #105 
		    #"met,mT2W,HTfrac,dPhi_JetMet,lepton_pT,dR_LepB,Chi2SNT,jet1_pT"  #106
		    #"met,lepton_pT,mlb_hemi,HTfrac,dPhi_JetMet,b1_pt" # 107


		    #t2bw	
                    #"met,mT2W,HTfrac,b1_pt,dPhi_JetMet,dR_LepB" #200
                    #"met,lepton_pT,mT2W,HTfrac,b1_pt,dPhi_JetMet,dR_LepB" #201
		    #t2tt	
                    #"met,mT2W,HTfrac,dPhi_JetMet,Chi2SNT" #202
                    #"met,mT2W,HTfrac,dPhi_JetMet,b1_pt" #203
                  ]



    NN_book_options_list = []
    NN_book_options_list.append("!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning")

    
    counter = 105 
    train_test = "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=EqualNumEvents"


    for NN_vars in NN_vars_list:
        for NN_book_options in NN_book_options_list:
            configfilename = getconfigname(counter,basedir)
            print configfilename
            backgroundfilelist = os.getcwd() + "/filelists/filelist_TTbar_skim.txt"
            signalsample = os.getcwd() + "/filelists/filelist_signal_skim.txt"
            makeconfigfile(configfilename,NN_book_options,NN_vars,train_test,backgroundfilelist,signalsample)
            counter = counter + 1
            


    
