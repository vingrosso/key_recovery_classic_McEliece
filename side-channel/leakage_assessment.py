import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import f_oneway
import argparse

def leakage_assessment(anova_test_output, traces_input, intermediate_values_input, mceliece, decapsulation_traces=True, lda=False):
    traces_matrix = np.load(traces_input)
    if mceliece == '8192128':
        if lda : 
            traces_matrix = traces_matrix["inv_traces"]
        else : traces_matrix = traces_matrix["traces_matrix"][:,:,:31350]
    else:
        traces_matrix = traces_matrix["traces_matrix"]
    print(traces_matrix.shape)

    hamming_weight_intermediate_values_matrix = np.load(intermediate_values_input)
    if lda : 
        hamming_weight_intermediate_values_matrix = hamming_weight_intermediate_values_matrix['inv_intermediate_values']
        hamming_weight_intermediate_values_matrix = np.array([[bin(intermediate_value).count('1') for intermediate_value in hamming_weight_intermediate_values_matrix[index_trace]] for index_trace in range(hamming_weight_intermediate_values_matrix.shape[0])])
    else:
        if mceliece == '8192128':
            hamming_weight_intermediate_values_matrix = hamming_weight_intermediate_values_matrix['intermediate_values_matrix'][:,:,:5]
        else:
            hamming_weight_intermediate_values_matrix = hamming_weight_intermediate_values_matrix['intermediate_values_matrix']
        hamming_weight_intermediate_values_matrix = np.array([[[bin(intermediate_value).count('1') for intermediate_value in hamming_weight_intermediate_values_alpha] for hamming_weight_intermediate_values_alpha in hamming_weight_intermediate_values_matrix[index_trace]] for index_trace in range(traces_matrix.shape[0])])
    print(hamming_weight_intermediate_values_matrix.shape)

    if not decapsulation_traces:
        if mceliece == '348864':
            nb_anova_test_traces = 20000
        elif mceliece == '8192128':
            if lda : nb_anova_test_traces = -1
            else : nb_anova_test_traces = 30000
        if lda : hamming_weight_intermediate_values_matrix =  hamming_weight_intermediate_values_matrix.reshape(hamming_weight_intermediate_values_matrix.shape[0]*hamming_weight_intermediate_values_matrix.shape[1])
        else : hamming_weight_intermediate_values_matrix =  hamming_weight_intermediate_values_matrix.reshape(hamming_weight_intermediate_values_matrix.shape[0]*hamming_weight_intermediate_values_matrix.shape[1],hamming_weight_intermediate_values_matrix.shape[-1])
        print(hamming_weight_intermediate_values_matrix.shape)
        if not lda : hamming_weight_intermediate_values_matrix = hamming_weight_intermediate_values_matrix[:nb_anova_test_traces,:]
        traces_matrix = traces_matrix.reshape(traces_matrix.shape[0]*traces_matrix.shape[1],traces_matrix.shape[-1])
        if not lda : traces_matrix = traces_matrix[:nb_anova_test_traces,:]
        print(traces_matrix.shape)

    # Compute the f-test for several intermediate values

    if lda : nb_val_int_to_test = 1
    elif not decapsulation_traces : nb_val_int_to_test = hamming_weight_intermediate_values_matrix.shape[1]
    else : 
        nb_val_int_to_test = hamming_weight_intermediate_values_matrix.shape[2]
    
    print(nb_val_int_to_test)

    nb_traces = traces_matrix.shape[0]
    if decapsulation_traces : 
        nb_samples = traces_matrix.shape[1]
        nb_alpha = hamming_weight_intermediate_values_matrix.shape[1]
    else : 
        nb_samples = traces_matrix.shape[1]
        nb_alpha = 1

    if decapsulation_traces : l_F_e_inv = np.zeros((nb_alpha, nb_val_int_to_test, nb_samples),dtype=float)
    else : l_F_e_inv = np.zeros((nb_val_int_to_test, nb_samples),dtype=float)

    print(l_F_e_inv.shape)

    for index_alpha in range(nb_alpha):
        for val_int_to_test in range(nb_val_int_to_test):
            nb_classes_max = np.max(hamming_weight_intermediate_values_matrix) + 1
            print(nb_classes_max)
            l_nb_traces_per_class = [0]*nb_classes_max

            # Get the number of traces in each class

            for index_trace in range(nb_traces):
                if lda : index_class = hamming_weight_intermediate_values_matrix[index_trace]
                elif decapsulation_traces : index_class = hamming_weight_intermediate_values_matrix[index_trace][index_alpha][val_int_to_test]
                else : index_class = hamming_weight_intermediate_values_matrix[index_trace][val_int_to_test]
                l_nb_traces_per_class[index_class] +=1

            assert sum(l_nb_traces_per_class) == nb_traces
            print(l_nb_traces_per_class)
            
            # Classes where the number of traces is equal to zero are removed from the test

            class_to_remove = []

            for index_list in range(len(l_nb_traces_per_class)):
                if not l_nb_traces_per_class[index_list]:
                    print('The class '+str(index_list)+' has been removed')
                    class_to_remove.append(index_list)

            class_to_remove.reverse()
            
            for index_to_remove in class_to_remove:
                del l_nb_traces_per_class[index_to_remove]

            # Create the arrays for each class stored in a list

            HW_val_int_classes = []

            for nb_traces_per_class in l_nb_traces_per_class:
                HW_val_int_classes.append(np.zeros((nb_traces_per_class, nb_samples), dtype="float"))

            # Traces classification

            nb_row = [0]*len(l_nb_traces_per_class)
            for index_trace in range(nb_traces):
                if lda : index_class = hamming_weight_intermediate_values_matrix[index_trace]
                elif decapsulation_traces : index_class = hamming_weight_intermediate_values_matrix[index_trace][index_alpha][val_int_to_test]
                else : index_class = hamming_weight_intermediate_values_matrix[index_trace][val_int_to_test]
                for index_to_remove in class_to_remove:
                    if index_class > index_to_remove :
                        index_class -= 1
                HW_val_int_classes[index_class][nb_row[index_class]][:] = traces_matrix[index_trace]
                
                nb_row[index_class]+=1

            #F-test

            F, p = f_oneway(*HW_val_int_classes)
            if lda : l_F_e_inv[val_int_to_test][:] = F
            elif decapsulation_traces : l_F_e_inv[index_alpha][val_int_to_test][:] = F
            else : l_F_e_inv[val_int_to_test][:] = F

    print(l_F_e_inv.shape)
    print(l_F_e_inv[0])
    np.savez_compressed(anova_test_output, F_e_inv_matrix=l_F_e_inv)

def plot_leakage_assessment(traces_input, anova_test_output, anova_plot_output, mceliece, all_alphas=False, alpha_traces=False, lda= False):
    traces_matrix = np.load(traces_input)
    if mceliece == '8192128':
        if lda : traces_matrix = traces_matrix["inv_traces"]
        else : traces_matrix = traces_matrix["traces_matrix"][:,:,:31350]
    else:
        traces_matrix = traces_matrix["traces_matrix"]
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({"font.size": 8})
    plt.rcParams["font.family"] = "serif" 

    fig, axes_anova = plt.subplots(figsize=(5.25,2.06))

    
    xlabel = "Time samples"
    ylabel = "F statistic"
    ylabel_trace = "Power consumption"

    dpi = 300
    labelspacing = 0.4
    columnspacing = 0.4
    
    nb_col_legend = 5
    cm_f_e_inv = mpl.colormaps['YlOrRd']
    colors_f_e_inv = cm_f_e_inv(np.arange(240))
    color_F_e_inv = colors_f_e_inv[-100]

    print(cm_f_e_inv.N)
    print(type(cm_f_e_inv))
    print(type(colors_f_e_inv))
    axes_anova.set_xlabel(xlabel)
    
    if not alpha_traces and not lda:
        if all_alphas : 
            xlim = 32000
            axes_anova.set_xlim(0,xlim)
        else : axes_anova.set_xlim(0,14600)
    axes_anova.tick_params(axis='x')     
    axes_anova.set_ylabel(ylabel, color=color_F_e_inv)
    axes_anova.tick_params(axis='y', labelcolor=color_F_e_inv)        

    axes_trace_F_e_inv = axes_anova.twinx()
    axes_anova.set_zorder(axes_trace_F_e_inv.get_zorder()+1)
    axes_anova.patch.set_visible(False)
    color_trace = '0.85'
    axes_trace_F_e_inv.tick_params(axis='y', labelsize=0, size=0)
    for ylabel_i in axes_trace_F_e_inv.get_yticklabels():
        ylabel_i.set_visible(False)
        ylabel_i.set_fontsize(0.0)
    if alpha_traces or lda : axes_trace_F_e_inv.plot(traces_matrix[0][0], color=color_trace)
    else : axes_trace_F_e_inv.plot(traces_matrix[0], color=color_trace)

    with np.load(anova_test_output) as f_test_matrices:
        if alpha_traces or lda:
            f_e_inv_matrix = f_test_matrices['F_e_inv_matrix']
        else:
            if not all_alphas : f_e_inv_matrix = np.concatenate((f_test_matrices['F_e_inv_matrix'][0,:-1],f_test_matrices['F_e_inv_matrix'][1,:-31]))
            else : 
                f_e_inv_matrix = np.zeros((3, f_test_matrices['F_e_inv_matrix'][0,:-1].shape[0], f_test_matrices['F_e_inv_matrix'][0,:-1].shape[1]),dtype=float)
                f_e_inv_matrix[0] = f_test_matrices['F_e_inv_matrix'][0,:-1]
                f_e_inv_matrix[1] = f_test_matrices['F_e_inv_matrix'][1,:-1]
                f_e_inv_matrix[2] = f_test_matrices['F_e_inv_matrix'][2,:-1]
        print(f_e_inv_matrix.shape)

    
    if f_e_inv_matrix.shape[0] == 1 : index_color_f_e_inv = -100
    else : index_color_f_e_inv = 0
    
    if not all_alphas : 
        index = 0
        for f_e_inv in f_e_inv_matrix[0:6:2,:]:
            if index > 1 : axes_anova.plot(f_e_inv, color=colors_f_e_inv[index_color_f_e_inv], label=r"$\alpha_0^{"+str(index)+r"} (g^2(\alpha_0))^{-1}$")
            elif index : axes_anova.plot(f_e_inv, color=colors_f_e_inv[index_color_f_e_inv], label=r"$\alpha_0 (g^2(\alpha_0))^{-1}$")
            else : axes_anova.plot(f_e_inv, color=colors_f_e_inv[index_color_f_e_inv], label=r"$(g^2(\alpha_0))^{-1}$")
            index_color_f_e_inv+=len(colors_f_e_inv)//f_e_inv_matrix.shape[0]
            index += 1
    else:
        cm_f_e_inv_1 = mpl.colormaps['Blues']
        cm_f_e_inv_2 = mpl.colormaps['Greens']
        cm_f_e_inv_3 = mpl.colormaps['Oranges']
        colors_f_e_inv = []
        colors_f_e_inv.append(cm_f_e_inv_1(np.arange(cm_f_e_inv_1.N)))
        colors_f_e_inv.append(cm_f_e_inv_2(np.arange(cm_f_e_inv_1.N)))
        colors_f_e_inv.append(cm_f_e_inv_3(np.arange(cm_f_e_inv_1.N)))
        index_f_e_inv = 0
        for f_e_inv in f_e_inv_matrix:
            index_color_f_e_inv = 21
            for f_e in f_e_inv:
                axes_anova.plot(f_e, color=colors_f_e_inv[index_f_e_inv][index_color_f_e_inv])
                index_color_f_e_inv += len(colors_f_e_inv[index_f_e_inv])//f_e_inv.shape[0]
            index_f_e_inv += 1

    if not alpha_traces :
        box = axes_anova.get_position()
        axes_anova.set_position([box.x0, box.y0, box.width*0.9, box.height])

        if all_alphas :
            cb1ax = fig.add_axes([0.845, box.y0 + 2*((box.height-0.1)/3 + 0.05), 0.035, (box.height-0.1)/3])
            cb2ax = fig.add_axes([0.845, box.y0 + (box.height-0.1)/3 + 0.05, 0.035, (box.height-0.1)/3])
            cb3ax = fig.add_axes([0.845, box.y0, 0.035, (box.height-0.1)*0.33])
            cbar1 = fig.colorbar(mpl.cm.ScalarMappable(cmap='Oranges'), cax=cb1ax, boundaries=np.linspace(21/256,1, num=100), ticks=[35/256, 0.985])
            cbar1.ax.set_yticklabels([r"$(g^2(\alpha_2))^{-1}$", r"$\alpha_2^{39} (g^2(\alpha_2))^{-1}$"])
            cbar1.ax.tick_params(size=0)
            cbar2 = fig.colorbar(mpl.cm.ScalarMappable(cmap='Greens'), cax=cb2ax, boundaries=np.linspace(21/256,1, num=100), ticks=[35/256, 0.985])
            cbar2.ax.set_yticklabels([r"$(g^2(\alpha_1))^{-1}$", r"$\alpha_1^{39} (g^2(\alpha_1))^{-1}$"])
            cbar2.ax.tick_params(size=0)
            cbar3 = fig.colorbar(mpl.cm.ScalarMappable(cmap='Blues'), cax=cb3ax, boundaries=np.linspace(21/256,1, num=100), ticks=[35/256, 0.985])
            cbar3.ax.set_yticklabels([r"$(g^2(\alpha_0))^{-1}$", r"$\alpha_0^{39} (g^2(\alpha_0))^{-1}$"])
            cbar3.ax.tick_params(size=0)
        else:
            cbar=fig.colorbar(mpl.cm.ScalarMappable(cmap='YlOrRd'), ax=axes_anova, ticks=[0,0.765,0.815,1])
            cbar.ax.tick_params(size=0)
            cbar.ax.set_yticklabels([r"$(g^2(\alpha_0))^{-1}$", r"$\alpha_0^{39} (g^2(\alpha_0))^{-1}$", r"$(g^2(\alpha_1))^{-1}$", r"$\alpha_1^{9} (g^2(\alpha_1))^{-1}$"])   
    plt.savefig(anova_plot_output, bbox_inches='tight',dpi=dpi)

if __name__ == '__main__': 
    parser = argparse.ArgumentParser()

    parser.add_argument('--lda', action='store_true', dest='lda', default=False)

    parser.add_argument('--me','--mceliece', action='store', dest='mceliece', required=True,
                        choices=['51220', '348864', '8192128'],help='McEliece security paramaters')
        
    leakage = parser.add_mutually_exclusive_group(required=True)
    leakage.add_argument('-l','--leakage_assessment', dest='leakage_assessment', action='store_true')
    leakage.add_argument('--nl','--no-leakage_assessment', dest='leakage_assessment', action='store_false')

    plot = parser.add_mutually_exclusive_group(required=True)
    plot.add_argument('-p','--plot', dest='plot', action='store_true')
    plot.add_argument('--np','--no-plot', dest='plot', action='store_false')

    args = parser.parse_args()
    
    if args.mceliece == '51220':
        traces_input = 'results/matrices/'+str(args.mceliece)+'/leakage_assessment/decapsulation_traces_matrix_'+str(args.mceliece)+'.npz'
        intermediate_values_input = 'results/matrices/'+str(args.mceliece)+'/leakage_assessment/intermediate_values_'+str(args.mceliece)+'.npz'
    if args.mceliece == '348864':
        traces_input = 'results/matrices/'+str(args.mceliece)+'/attack/alpha_traces_'+str(args.mceliece)+'.npz'
        intermediate_values_input = 'results/matrices/'+str(args.mceliece)+'/attack/intermediate_values_'+str(args.mceliece)+'.npz'
    if args.mceliece == '8192128':
        traces_input = 'results/matrices/'+str(args.mceliece)+'/attack/alpha_traces_'+str(args.mceliece)+'_0.npz'
        intermediate_values_input = 'results/matrices/'+str(args.mceliece)+'/attack/intermediate_values_'+str(args.mceliece)+'_0.npz'

    anova_test_output = 'results/statistical_tests/'+str(args.mceliece)+'/anova_test_'+str(args.mceliece)+'.npz'
    anova_plot_output = 'results/statistical_tests/'+str(args.mceliece)+'/anova_plot_'+str(args.mceliece)+'.pdf'

    if args.mceliece == '51220':
        if args.leakage_assessment : leakage_assessment(anova_test_output, traces_input, intermediate_values_input, args.mceliece)
        if args.plot : plot_leakage_assessment(traces_input, anova_test_output, anova_plot_output, args.mceliece, True, False)
    if args.mceliece == '348864':
        if args.leakage_assessment : leakage_assessment(anova_test_output, traces_input, intermediate_values_input, args.mceliece, False)
        if args.plot : plot_leakage_assessment(traces_input, anova_test_output, anova_plot_output, args.mceliece, False, True)
    if args.mceliece == '8192128':
        if args.lda :
            intermediate_values_input = 'results/matrices/'+str(args.mceliece)+'/attack/lda_training_coefficient_intermediate_values_'+str(args.mceliece)+'.npz'
            traces_input = 'results/matrices/'+str(args.mceliece)+'/attack/lda_training_traces_'+str(args.mceliece)+'.npz'
            anova_test_output = 'results/statistical_tests/'+str(args.mceliece)+'/lda_anova_test_'+str(args.mceliece)+'.npz'
            anova_plot_output = 'results/statistical_tests/'+str(args.mceliece)+'/lda_anova_plot_'+str(args.mceliece)+'.pdf'
            if args.leakage_assessment : leakage_assessment(anova_test_output, traces_input, intermediate_values_input, args.mceliece, False, args.lda)
            if args.plot : plot_leakage_assessment(traces_input, anova_test_output, anova_plot_output, args.mceliece, False, False, args.lda)
        if args.leakage_assessment and not args.lda : leakage_assessment(anova_test_output, traces_input, intermediate_values_input, args.mceliece, False)
        if args.plot and not args.lda : plot_leakage_assessment(traces_input, anova_test_output, anova_plot_output, args.mceliece, False, True)