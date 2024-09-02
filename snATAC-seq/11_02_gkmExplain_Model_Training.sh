# Example of training of gkmExplain Models
lsgkm-master/src/gkmtrain -T 16 FA_file/EX_Human_Fold0_Pos_Train.fa FA_file/EX_Human_Fold0_Neg_Train.fa gkmExplain_Model/EX_Human_Fold0
lsgkm-master/src/gkmtrain -T 16 FA_file/EX_Human_Fold1_Pos_Train.fa FA_file/EX_Human_Fold1_Neg_Train.fa gkmExplain_Model/EX_Human_Fold1
lsgkm-master/src/gkmtrain -T 16 FA_file/EX_Human_Fold2_Pos_Train.fa FA_file/EX_Human_Fold2_Neg_Train.fa gkmExplain_Model/EX_Human_Fold2
lsgkm-master/src/gkmtrain -T 16 FA_file/EX_Human_Fold3_Pos_Train.fa FA_file/EX_Human_Fold3_Neg_Train.fa gkmExplain_Model/EX_Human_Fold3
lsgkm-master/src/gkmtrain -T 16 FA_file/EX_Human_Fold4_Pos_Train.fa FA_file/EX_Human_Fold4_Neg_Train.fa gkmExplain_Model/EX_Human_Fold4
lsgkm-master/src/gkmtrain -T 16 FA_file/EX_Human_Fold5_Pos_Train.fa FA_file/EX_Human_Fold5_Neg_Train.fa gkmExplain_Model/EX_Human_Fold5
lsgkm-master/src/gkmtrain -T 16 FA_file/EX_Human_Fold6_Pos_Train.fa FA_file/EX_Human_Fold6_Neg_Train.fa gkmExplain_Model/EX_Human_Fold6
lsgkm-master/src/gkmtrain -T 16 FA_file/EX_Human_Fold7_Pos_Train.fa FA_file/EX_Human_Fold7_Neg_Train.fa gkmExplain_Model/EX_Human_Fold7
lsgkm-master/src/gkmtrain -T 16 FA_file/EX_Human_Fold8_Pos_Train.fa FA_file/EX_Human_Fold8_Neg_Train.fa gkmExplain_Model/EX_Human_Fold8
lsgkm-master/src/gkmtrain -T 16 FA_file/EX_Human_Fold9_Pos_Train.fa FA_file/EX_Human_Fold9_Neg_Train.fa gkmExplain_Model/EX_Human_Fold9





# Example of predict HAR scores using gkmExplain Models (Reference Genome)
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Ref.fa ../gkmExplain_Model/EX_Human_Fold0.model.txt Res/EX_Human_Fold0_Ref.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Ref.fa ../gkmExplain_Model/EX_Human_Fold1.model.txt Res/EX_Human_Fold1_Ref.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Ref.fa ../gkmExplain_Model/EX_Human_Fold2.model.txt Res/EX_Human_Fold2_Ref.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Ref.fa ../gkmExplain_Model/EX_Human_Fold3.model.txt Res/EX_Human_Fold3_Ref.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Ref.fa ../gkmExplain_Model/EX_Human_Fold4.model.txt Res/EX_Human_Fold4_Ref.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Ref.fa ../gkmExplain_Model/EX_Human_Fold5.model.txt Res/EX_Human_Fold5_Ref.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Ref.fa ../gkmExplain_Model/EX_Human_Fold6.model.txt Res/EX_Human_Fold6_Ref.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Ref.fa ../gkmExplain_Model/EX_Human_Fold7.model.txt Res/EX_Human_Fold7_Ref.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Ref.fa ../gkmExplain_Model/EX_Human_Fold8.model.txt Res/EX_Human_Fold8_Ref.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Ref.fa ../gkmExplain_Model/EX_Human_Fold9.model.txt Res/EX_Human_Fold9_Ref.txt





# Example of predict HAR scores using gkmExplain Models (SNC Altered Genome)
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Anc.fa ../gkmExplain_Model/EX_Human_Fold0.model.txt Res/EX_Human_Fold0_Anc.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Anc.fa ../gkmExplain_Model/EX_Human_Fold1.model.txt Res/EX_Human_Fold1_Anc.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Anc.fa ../gkmExplain_Model/EX_Human_Fold2.model.txt Res/EX_Human_Fold2_Anc.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Anc.fa ../gkmExplain_Model/EX_Human_Fold3.model.txt Res/EX_Human_Fold3_Anc.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Anc.fa ../gkmExplain_Model/EX_Human_Fold4.model.txt Res/EX_Human_Fold4_Anc.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Anc.fa ../gkmExplain_Model/EX_Human_Fold5.model.txt Res/EX_Human_Fold5_Anc.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Anc.fa ../gkmExplain_Model/EX_Human_Fold6.model.txt Res/EX_Human_Fold6_Anc.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Anc.fa ../gkmExplain_Model/EX_Human_Fold7.model.txt Res/EX_Human_Fold7_Anc.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Anc.fa ../gkmExplain_Model/EX_Human_Fold8.model.txt Res/EX_Human_Fold8_Anc.txt
../lsgkm-master/src/gkmpredict -T 16 HAR_fa_file/HAR_Anc.fa ../gkmExplain_Model/EX_Human_Fold9.model.txt Res/EX_Human_Fold9_Anc.txt
