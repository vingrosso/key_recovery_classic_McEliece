import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal
from math import floor
from math import ceil
import argparse
from tqdm import tqdm 

def profiling_phase(LDA_traces, coefficient_intermediate_values, array_templates, maxima_minima, multivariate):
    
    inv_traces = np.load(LDA_traces)
    loop_traces = inv_traces['loop_traces']
    inv_traces = inv_traces['inv_traces']

    list_traces = [inv_traces, loop_traces]

    inv_intermediate_values = np.load(coefficient_intermediate_values)
    loop_intermediate_values = inv_intermediate_values['loop_intermediate_values']
    inv_intermediate_values = inv_intermediate_values['inv_intermediate_values']
    print(inv_intermediate_values.shape)
    print(loop_intermediate_values.shape)
    inv_hw_intermediate_values = np.array([[bin(intermediate_value).count('1') for intermediate_value in inv_intermediate_values[index_trace]] for index_trace in range(inv_intermediate_values.shape[0])])
    print(inv_hw_intermediate_values.shape)
    loop_hw_intermediate_values = np.array([[[bin(intermediate_value).count('1') for intermediate_value in hamming_weight_intermediate_values_alpha] for hamming_weight_intermediate_values_alpha in loop_intermediate_values[index_trace]] for index_trace in range(loop_intermediate_values.shape[0])])
    print(loop_hw_intermediate_values.shape)

    list_Hamming_weight = [inv_hw_intermediate_values, loop_hw_intermediate_values]  
    
    list_templates = np.zeros((2,),dtype=np.ndarray)
    max_values = np.zeros((2),dtype=float)
    min_values = np.zeros((2),dtype=float)

    index_templates = 0
    for traces_matrix, hamming_weight_intermediate_values_matrix in zip(list_traces, list_Hamming_weight):
        assert traces_matrix.shape[0] == hamming_weight_intermediate_values_matrix.shape[0]
        nb_decapsulation_traces = traces_matrix.shape[0]
        nb_alpha_per_traces = traces_matrix.shape[1]
        if len(traces_matrix.shape) == 4 : nb_loop_per_alpha = traces_matrix.shape[2]
        nb_samples = traces_matrix.shape[-1]

        nb_classes_max = np.max(hamming_weight_intermediate_values_matrix) + 1
        print(nb_classes_max)
        l_nb_traces_per_class = [0]*nb_classes_max

        # Get the number of traces in each class

        for index_trace in range(nb_decapsulation_traces):
            for index_alpha in range(nb_alpha_per_traces):
                if len(traces_matrix.shape) == 4:
                    for index_loop in range(nb_loop_per_alpha):
                        index_class = hamming_weight_intermediate_values_matrix[index_trace][index_alpha][index_loop]
                        l_nb_traces_per_class[index_class] +=1
                else:
                    index_class = hamming_weight_intermediate_values_matrix[index_trace][index_alpha]
                    l_nb_traces_per_class[index_class] +=1

        if len(traces_matrix.shape) == 4 : assert sum(l_nb_traces_per_class) == nb_decapsulation_traces*nb_alpha_per_traces*nb_loop_per_alpha
        else : assert sum(l_nb_traces_per_class) == nb_decapsulation_traces*nb_alpha_per_traces

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
            HW_val_int_classes.append(np.zeros((nb_traces_per_class,nb_samples), dtype="float"))

        # Traces classification

        nb_row = [0]*len(l_nb_traces_per_class)
        if len(traces_matrix.shape) == 4:
            for index_trace in range(nb_decapsulation_traces):
                for index_alpha in range(nb_alpha_per_traces):
                    for index_loop in range(nb_loop_per_alpha):
                        index_class = hamming_weight_intermediate_values_matrix[index_trace][index_alpha][index_loop]
                        for index_to_remove in class_to_remove:
                            if index_class > index_to_remove :
                                index_class -= 1
                        HW_val_int_classes[index_class][nb_row[index_class]][:] = traces_matrix[index_trace][index_alpha][index_loop]
                        nb_row[index_class]+=1
        else:
            for index_trace in range(nb_decapsulation_traces):
                for index_alpha in range(nb_alpha_per_traces):
                    index_class = hamming_weight_intermediate_values_matrix[index_trace][index_alpha]
                    for index_to_remove in class_to_remove:
                        if index_class > index_to_remove :
                            index_class -= 1
                    HW_val_int_classes[index_class][nb_row[index_class]][:] = traces_matrix[index_trace][index_alpha]
                    nb_row[index_class]+=1

        # Compute the templates (univariate normal distribution)

        index_poi = 0
        templates = np.zeros((len(HW_val_int_classes)), dtype=type(multivariate_normal()))
        index_HW_class = 0
        for HW_class in HW_val_int_classes:
            if HW_class[:,index_poi].max() > max_values[index_templates] : max_values[index_templates] = HW_class[:,index_poi].max() + 5
            if HW_class[:,index_poi].min() < min_values[index_templates] : min_values[index_templates] = HW_class[:,index_poi].min() - 5
            if multivariate :
                templates[index_HW_class] = multivariate_normal(np.mean(HW_class[:,:], axis=0), np.cov(HW_class[:,:].T))
            else :
                templates[index_HW_class] = multivariate_normal(np.mean(HW_class[:,index_poi], axis=0), np.cov(HW_class[:,index_poi].T))
            index_HW_class += 1
        list_templates[index_templates] = templates
        index_templates += 1

    np.savez_compressed(array_templates, templates_inv=list_templates[0], templates_loop=list_templates[1])
    np.savez_compressed(maxima_minima, maxima=max_values, minima=min_values)

def plot_profiling_phase(templates_plot, array_templates, maxima_minima):
    
    list_templates = np.load(array_templates,allow_pickle=True)
    list_templates = list(list_templates.values())
    maxima = np.load(maxima_minima)
    minima = maxima['minima']
    maxima = maxima['maxima']

    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({"font.size": 8})
    plt.rcParams["font.family"] = "serif"

    fig, axes_templates = plt.subplots(nrows = len(list_templates))
    
    labelspacing = 0.4
    columnspacing = 0.4

    index_ax = 0
    for ax,templates,maximum,minimum in zip(axes_templates.flat, list_templates, maxima, minima):
        print(minimum, maximum)
        x = np.linspace(floor(minimum), ceil(maximum), 2000)
        ax.set_ylabel("Probability")
        ax.tick_params(axis='y')
        ax.set_xlabel("Power consumption")
        ax.tick_params(axis='x')
        if index_ax:
            for index_HW in range(len(templates)):
                ax.plot(x, templates[index_HW].pdf(x), label=f"wt={index_HW}")
        else:
            for index_HW in range(len(templates)):
                label_HW = index_HW + 1
                ax.plot(x, templates[index_HW].pdf(x), label=f"wt={label_HW}")
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        nb_col_legend = len(templates)//2 + 1
        handles, labels = ax.get_legend_handles_labels()
        order = []
        for index_column in range(nb_col_legend):
            if index_column >= len(labels)%nb_col_legend : nb_row = len(labels)//nb_col_legend
            else : nb_row = len(labels)//nb_col_legend + 1
            for index_row in range(nb_row):
                order.append(index_column + index_row * nb_col_legend)
        ax.legend(handles = [handles[i] for i in order], labels = [labels[i] for i in order],loc='lower center', bbox_to_anchor=(0.5, -0.35),fancybox=True, 
                  ncol=nb_col_legend, columnspacing = columnspacing, labelspacing = labelspacing)
        index_ax += 1
    
    pos = axes_templates[1].get_position().get_points()
    pos[0][1] -= 0.052
    pos[1][1] -= 0.052
    axes_templates[1].set_position(mpl.transforms.Bbox(pos))
    plt.show()

def matching_phase(array_probabilities, LDA_traces, array_templates, multivariate, mceliece):
    
    templates_inv_matrix = np.load(array_templates, allow_pickle=True)
    templates_loop_matrix = templates_inv_matrix['templates_loop']
    templates_inv_matrix = templates_inv_matrix['templates_inv']
    
    inv_test_traces_matrix = np.load(LDA_traces)
    loop_test_traces_matrix = inv_test_traces_matrix['loop_traces']
    inv_test_traces_matrix = inv_test_traces_matrix['inv_traces']
    
    print(inv_test_traces_matrix.shape)
    print(loop_test_traces_matrix.shape)

    nb_HW_classes = loop_test_traces_matrix.shape[3]

    nb_decapsulation_traces = inv_test_traces_matrix.shape[0]
    nb_alpha_per_trace = inv_test_traces_matrix.shape[1]
    nb_coefficients_per_column = loop_test_traces_matrix.shape[2] + 1

    nb_loop_traces = loop_test_traces_matrix.shape[2]

    probabilities_matrix = np.zeros((nb_decapsulation_traces, nb_alpha_per_trace, nb_coefficients_per_column, nb_HW_classes),dtype=float)
    log_probabilities_matrix = np.full((nb_decapsulation_traces, nb_alpha_per_trace,nb_coefficients_per_column, nb_HW_classes),fill_value=np.NINF,dtype=float)
    
    print(probabilities_matrix.shape)

    for index_trace in tqdm(range(nb_decapsulation_traces)):
        for index_alpha in range(nb_alpha_per_trace):
            for index_HW_class in range(nb_HW_classes-1):
                if multivariate:
                    probabilities_matrix[index_trace][index_alpha][0][index_HW_class+1] = templates_inv_matrix[index_HW_class].pdf(inv_test_traces_matrix[index_trace,index_alpha,:])
                    log_probabilities_matrix[index_trace][index_alpha][0][index_HW_class+1] = templates_inv_matrix[index_HW_class].logpdf(inv_test_traces_matrix[index_trace,index_alpha,:])
                else:
                    probabilities_matrix[index_trace][index_alpha][0][index_HW_class+1] = templates_inv_matrix[index_HW_class].pdf(inv_test_traces_matrix[index_trace][index_alpha][0])
                    log_probabilities_matrix[index_trace][index_alpha][0][index_HW_class+1] = templates_inv_matrix[index_HW_class].logpdf(inv_test_traces_matrix[index_trace][index_alpha][0])
            if mceliece != '8192128':
                for index_loop_trace in range(nb_loop_traces):
                    for index_HW_class in range(nb_HW_classes):
                        if multivariate:
                            probabilities_matrix[index_trace][index_alpha][index_loop_trace+1][index_HW_class] = templates_loop_matrix[index_HW_class].pdf(loop_test_traces_matrix[index_trace, index_alpha, index_loop_trace, :])
                            log_probabilities_matrix[index_trace][index_alpha][index_loop_trace+1][index_HW_class] = templates_loop_matrix[index_HW_class].logpdf(loop_test_traces_matrix[index_trace, index_alpha, index_loop_trace, :]) 
                        else:
                            probabilities_matrix[index_trace][index_alpha][index_loop_trace+1][index_HW_class] = templates_loop_matrix[index_HW_class].pdf(loop_test_traces_matrix[index_trace][index_alpha][index_loop_trace][0])
                            log_probabilities_matrix[index_trace][index_alpha][index_loop_trace+1][index_HW_class] = templates_loop_matrix[index_HW_class].logpdf(loop_test_traces_matrix[index_trace][index_alpha][index_loop_trace][0])
    if mceliece == '8192128':
        for index_trace in range(loop_test_traces_matrix.shape[0]):
            for index_alpha in range(loop_test_traces_matrix.shape[1]):
                for index_loop_trace in range(loop_test_traces_matrix.shape[2]):
                    for index_HW_class in range(nb_HW_classes):
                        probabilities_matrix[index_trace][index_alpha][index_loop_trace+1][index_HW_class] = templates_loop_matrix[index_HW_class].pdf(loop_test_traces_matrix[index_trace][index_alpha][index_loop_trace][0])
                        log_probabilities_matrix[index_trace][index_alpha][index_loop_trace+1][index_HW_class] = templates_loop_matrix[index_HW_class].logpdf(loop_test_traces_matrix[index_trace][index_alpha][index_loop_trace][0])

    np.savez_compressed(array_probabilities, probabilities=probabilities_matrix, log_probabilities=log_probabilities_matrix)

def check_results(array_probabilities, coefficient_intermediate_values, multivariate, mceliece):
    probabilities_matrix = np.load(array_probabilities)
    probabilities_matrix = probabilities_matrix['probabilities']
        
    inv_intermediate_values = coefficient_intermediate_values[0]
    loop_intermediate_values = coefficient_intermediate_values[1]
    print(inv_intermediate_values.shape)
    print(loop_intermediate_values.shape)
    inv_hw_intermediate_values = np.array([[bin(intermediate_value).count('1') for intermediate_value in inv_intermediate_values[index_trace]] for index_trace in range(inv_intermediate_values.shape[0])])
    print(inv_hw_intermediate_values.shape)
    loop_hw_intermediate_values = np.array([[[bin(intermediate_value).count('1') for intermediate_value in hamming_weight_intermediate_values_alpha] for hamming_weight_intermediate_values_alpha in loop_intermediate_values[index_trace]] for index_trace in range(loop_intermediate_values.shape[0])])
    print(loop_hw_intermediate_values.shape)

    print(probabilities_matrix.shape)
    print(np.argmax(probabilities_matrix, axis=-1).shape)

    nb_traces_test_set = len(probabilities_matrix)

    nb_decapsulation_traces = inv_intermediate_values.shape[0]
    nb_alpha_per_trace = inv_intermediate_values.shape[1]
    nb_loop_traces = loop_intermediate_values.shape[2]

    HW_test = np.argmax(probabilities_matrix, axis=-1)

    counter_inv = 0
    for index_trace in range(nb_decapsulation_traces):
        for index_alpha in range(nb_alpha_per_trace): 
            if inv_hw_intermediate_values[index_trace][index_alpha] == HW_test[index_trace][index_alpha][0]: 
                counter_inv += 1
    print(nb_decapsulation_traces*nb_alpha_per_trace)
    print(counter_inv)
    print(nb_decapsulation_traces*nb_alpha_per_trace - counter_inv)
    
    if mceliece == '8192128':
        counter_loop = 0
        for index_trace in range(loop_intermediate_values.shape[0]):
            for index_alpha in range(nb_alpha_per_trace): 
                for index_loop in range(nb_loop_traces):
                    if loop_hw_intermediate_values[index_trace][index_alpha][index_loop] == HW_test[index_trace][index_alpha][index_loop+1]:
                        counter_loop += 1
    else :
        counter_loop = 0
        for index_trace in range(nb_decapsulation_traces):
            for index_alpha in range(nb_alpha_per_trace): 
                for index_loop in range(nb_loop_traces):
                    if loop_hw_intermediate_values[index_trace][index_alpha][index_loop] == HW_test[index_trace][index_alpha][index_loop+1]: 
                        counter_loop += 1
    print(loop_intermediate_values.shape[0]*nb_alpha_per_trace*nb_loop_traces)
    print(counter_loop)
    print(loop_intermediate_values.shape[0]*nb_alpha_per_trace*nb_loop_traces - counter_loop)

if __name__ == "__main__":
    import argparse
 
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-p','--profiling', action='store_true', dest='profiling')
    
    parser.add_argument('-d','--plot', action='store_true', dest='plot')
    
    parser.add_argument('-m','--matching', action='store_true', dest='matching')

    parser.add_argument('-c','--check', action='store_true', dest='check_results')
    
    parser.add_argument('--ml','--multivariate', action='store_true', dest='multivariate', default=False)
    
    parser.add_argument('--me','--mceliece', action='store', dest='mceliece', required=True,
                        choices=['51220', '348864', '8192128'],help='McEliece security paramaters')
    

    args = parser.parse_args()
    
    LDA_training_traces = 'results/matrices/'+str(args.mceliece)+'/attack/lda_training_traces_'+str(args.mceliece)+'.npz'
    lda_training_coefficient_intermediate_values = 'results/matrices/'+str(args.mceliece)+'/attack/lda_training_coefficient_intermediate_values_'+str(args.mceliece)+'.npz'

    if args.multivariate:
        array_templates = 'results/matrices/'+str(args.mceliece)+'/attack/multivariate_templates_'+str(args.mceliece)+'.npz'
        maxima_minima = 'results/matrices/'+str(args.mceliece)+'/attack/multivariate_maxima_minima_'+str(args.mceliece)+'.npz'
        array_probabilities ='results/matrices/'+str(args.mceliece)+'/attack/multivariate_probabilities_'+str(args.mceliece)+'.npz'
    else:
        array_templates = 'results/matrices/'+str(args.mceliece)+'/attack/templates_'+str(args.mceliece)+'.npz'
        maxima_minima = 'results/matrices/'+str(args.mceliece)+'/attack/maxima_minima_'+str(args.mceliece)+'.npz'
        array_probabilities ='results/matrices/'+str(args.mceliece)+'/attack/probabilities_'+str(args.mceliece)+'.npz'

    templates_plot = 'results/matrices/'+str(args.mceliece)+'/attack/templates_'+str(args.mceliece)+'.pdf'

    LDA_test_traces = 'results/matrices/'+str(args.mceliece)+'/attack/lda_test_traces_'+str(args.mceliece)+'.npz'
    LDA_test_coefficient_intermediates_values = 'results/matrices/'+str(args.mceliece)+'/attack/lda_test_coefficient_intermediates_values_'+str(args.mceliece)+'.npz'

    if args.check_results and args.mceliece != '8192128' :
        coefficient_intermediate_values = 'results/matrices/'+str(args.mceliece)+'/attack/coefficient_intermediates_values_'+str(args.mceliece)+'.npz'
        inv_intermediate_values = np.load(coefficient_intermediate_values)
        loop_intermediate_values = inv_intermediate_values['loop_intermediate_values']
        inv_intermediate_values = inv_intermediate_values['inv_intermediate_values']
    
    if args.mceliece == '51220':
        nb_training_traces_inv = 200
        nb_training_traces_loop = 200
        nb_test_traces_inv = nb_training_traces_inv
        nb_test_traces_loop = nb_training_traces_loop
    if args.mceliece == '348864':
        nb_training_traces_inv = 500
        nb_training_traces_loop = 300
        nb_test_traces_inv = 500
        nb_test_traces_loop = 500
    if args.mceliece == '8192128':
        nb_test_traces_inv = 700

    if args.profiling : profiling_phase(LDA_training_traces, lda_training_coefficient_intermediate_values, array_templates, maxima_minima, args.multivariate)
    if args.plot : plot_profiling_phase(templates_plot, array_templates, maxima_minima)
    if args.matching : matching_phase(array_probabilities, LDA_test_traces, array_templates, args.multivariate,  args.mceliece)
    if args.check_results : 
        if args.mceliece == '8192128' :
            list_coefficient_intermediate_values = list('results/matrices/'+str(args.mceliece)+'/attack/coefficient_intermediates_values_'+str(args.mceliece)+'_'+str(nb_file)+'.npz' for nb_file in range(5))
            inv_intermediate_values = np.load( list_coefficient_intermediate_values[0])
            loop_intermediate_values = inv_intermediate_values['loop_intermediate_values'][49:51]
            inv_intermediate_values = inv_intermediate_values['inv_intermediate_values']
            for path_to_matrix in list_coefficient_intermediate_values[1:] :
                with np.load(path_to_matrix) as coefficient_intermediate_values :
                    inv_intermediate_values = np.concatenate((inv_intermediate_values, coefficient_intermediate_values['inv_intermediate_values']), axis=0)
            print(inv_intermediate_values.shape)
            print(loop_intermediate_values.shape)
            check_results(array_probabilities, (inv_intermediate_values[nb_test_traces_inv:], loop_intermediate_values), args.multivariate, args.mceliece)
        else :
            check_results(array_probabilities, (inv_intermediate_values[nb_test_traces_inv:], loop_intermediate_values[nb_test_traces_loop:]), args.multivariate, args.mceliece)