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
                         		  "m3b,mlb_hemi,njets,dR_LepB,dPhi_JetMet,met,mT2W,lepton_pT,b1_pt" # t2bw025-50
                                        #    "m3b,mlb_hemi,njets,dPhi_JetMet,HT,met,mT2W,lepton_pT,jet1_pT", #t2bw075
		                        #     "met,HTfrac,njets,dPhi_JetMet,lepton_pT,mlb_hemi,b1_pt", #t2tt ofshell
                                        #      "met,mT2W,HTfrac,njets,dPhi_JetMet,lepton_pT,dR_LepB,Chi2SNT,jet1_pT", # t2tt onshell
		   ]



    NN_book_options_list = []
    NN_book_options_list.append("!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning")

    
    counter = 0 
    train_test = "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=EqualNumEvents"


    for NN_vars in NN_vars_list:
        for NN_book_options in NN_book_options_list:
            configfilename = getconfigname(counter,basedir)
            print configfilename
            backgroundfilelist = os.getcwd() + "/filelists/filelist_TTbar_skim.txt"
            signalsample = os.getcwd() + "/filelists/filelist_signal_skim.txt"
            makeconfigfile(configfilename,NN_book_options,NN_vars,train_test,backgroundfilelist,signalsample)
            counter = counter + 1
            
