hdist = fread('evaluation/iteration2_bubble_phase_evaluation/haplo_edit_dist_k63_a3_l23_kminimap15_w1_f0.1_z500_HG00733.data')
hdist.filt <- hdist[dist_h1!='None' & dist_h2!='None' & pred_type!='None']

ggplot(hdist.filt, aes(x=dist_h1, y=dist_h2, color=pred_type))+geom_point()
