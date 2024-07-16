import numpy as np
import galois as gf
from scipy.signal import correlate
from tqdm import tqdm

def decapsulation_traces_slicing(alpha_traces_output, list_decapsulation_traces_input, mceliece):
    if mceliece == '51220':
        offset = 0
        end_trace = -1
    if mceliece == '348864':
        offset = 15680
        end_trace = 4_192_000
    if mceliece == '8192128':
        offset = 0
        end_trace = -1
    nb_decapsulation_traces = 0

    for path_to_matrix in list_decapsulation_traces_input:
        with np.load(path_to_matrix) as decapsulation_traces_matrix:
            decapsulation_traces_matrix = decapsulation_traces_matrix["traces_matrix"][:,offset:end_trace]
            nb_decapsulation_traces += decapsulation_traces_matrix.shape[0]
            decapsulation_trace = decapsulation_traces_matrix[0]

    print(nb_decapsulation_traces)
    print(decapsulation_traces_matrix.shape)
    decapsulation_traces_matrix = np.load(list_decapsulation_traces_input[0])
    decapsulation_traces_matrix = decapsulation_traces_matrix["traces_matrix"][:,offset:end_trace]
    print(decapsulation_traces_matrix.shape[0])
    decapsulation_trace = decapsulation_traces_matrix[0]
    autocorrelation = correlate(decapsulation_trace,decapsulation_trace)
    period_trace = np.argmax(autocorrelation[(autocorrelation.size//2)+1:])+1 
    print(period_trace)
    nb_alpha_per_trace = decapsulation_trace.size//period_trace
    print(nb_alpha_per_trace)

    if mceliece == '8192128' :
        alpha_traces_matrix = np.zeros((nb_decapsulation_traces, nb_alpha_per_trace, period_trace-29_310), dtype="float")
    else : alpha_traces_matrix = np.zeros((nb_decapsulation_traces, nb_alpha_per_trace, period_trace),dtype="float")

    index_row = 0
    for path_to_matrix in list_decapsulation_traces_input:
        with np.load(path_to_matrix) as decapsulation_traces_matrix:
            decapsulation_traces_matrix = decapsulation_traces_matrix["traces_matrix"][:,offset:end_trace]
            for index_decapsulation_trace in tqdm(range(decapsulation_traces_matrix.shape[0])):
                for index_alpha_trace in range(nb_alpha_per_trace):
                    if mceliece == '8192128' :
                        alpha_traces_matrix[index_row][index_alpha_trace][:] = decapsulation_traces_matrix[index_decapsulation_trace][index_alpha_trace*period_trace:(index_alpha_trace+1)*period_trace-29_310]
                    else : alpha_traces_matrix[index_row][index_alpha_trace][:] = decapsulation_traces_matrix[index_decapsulation_trace][index_alpha_trace*period_trace:(index_alpha_trace+1)*period_trace]
                index_row += 1
    if mceliece == '8192128' : np.savez_compressed(alpha_traces_output, traces_matrix=alpha_traces_matrix);
    else : np.savez_compressed(alpha_traces_output, traces_matrix=alpha_traces_matrix)

def intermediate_values_computation(intermediates_values_output, alpha_traces_input, parameters_input, mceliece):
    parameters_matrices = np.load(parameters_input)

    g_matrix = parameters_matrices["g_matrix"]
    L_matrix = parameters_matrices["L_matrix"]
    ct_matrix = parameters_matrices["ct_matrix"]
    
    if mceliece == "51220":
        m=9
        t=20
        f_x = "x^9 + x^1 + 1"
        beginning_alpha = 0
    elif mceliece == "348864":
        m=12
        t=64
        f_x = "x^12 + x^3 + 1"
        beginning_alpha = 1
    elif mceliece == "8192128":
        m=13
        t=128
        f_x = "x^13 + x^4 + x^3 + x + 1"
        beginning_alpha = 0

    K = gf.GF(2**m, repr="poly", irreducible_poly=f_x)

    nb_val_int_before_loop = 1
    nb_val_int_loop = 1
    nb_iterations = 2*t
    nb_columns = nb_val_int_before_loop + nb_iterations*nb_val_int_loop

    
    if mceliece == "8192128":
        beginning_decapsulation = 400
        nb_decapsulation_traces = 200
        nb_alpha_traces = 74
    else :
        alpha_traces_matrix = np.load(alpha_traces_input)
        alpha_traces_matrix = alpha_traces_matrix['traces_matrix']
        nb_decapsulation_traces = alpha_traces_matrix.shape[0]
        nb_alpha_traces = alpha_traces_matrix.shape[1]
    
    intermediate_values_matrix = np.zeros((nb_decapsulation_traces, nb_alpha_traces, nb_columns),dtype=np.uint16)
    print(intermediate_values_matrix.shape)
    print(len(range(beginning_alpha, nb_alpha_traces+beginning_alpha)))
    print(beginning_alpha)
    print(nb_alpha_traces+beginning_alpha)
    for index_decapsulation_trace in tqdm(range(beginning_decapsulation, beginning_decapsulation+nb_decapsulation_traces)):
        for index_alpha_trace in range(beginning_alpha, nb_alpha_traces+beginning_alpha):
            g = gf.Poly(np.concatenate(([1], ((np.uint16(g_matrix[index_decapsulation_trace][-3::-2])<< 8)+g_matrix[index_decapsulation_trace][-4::-2]))), field=K)
            alpha = K([L_matrix[index_decapsulation_trace][2*index_alpha_trace] + (L_matrix[index_decapsulation_trace][2*index_alpha_trace+1] << 8)])

            e=g(alpha)
            e_inv = K([1])/(g(alpha)**2)
            intermediate_values_matrix[index_decapsulation_trace-beginning_decapsulation][index_alpha_trace-beginning_alpha][0] = e_inv

            for iterations in range(nb_iterations):
                e_inv = e_inv*alpha
                intermediate_values_matrix[index_decapsulation_trace-beginning_decapsulation][index_alpha_trace-beginning_alpha][iterations+1] = e_inv
    np.savez_compressed(intermediates_values_output, intermediate_values_matrix=intermediate_values_matrix)

def alpha_traces_slicing(coefficient_traces, coefficient_intermediate_values, anova_matrix, alpha_traces, intermediate_values, mceliece):
    if mceliece == '51220':
        m=9
        t=20
        n=512
    if mceliece == '348864':
        m=12
        t=64
        n=3488
    if mceliece == '8192128':
        m=13
        t=128
        n=8192
    anova_matrix = np.load(anova_matrix)
    if mceliece == '8192128':
        anova_matrix = anova_matrix['F_e_inv_matrix'][:3,:8000]
        print(anova_matrix.shape)
    else :
        anova_matrix = anova_matrix['F_e_inv_matrix']
    traces_separation_list = []
    index_anova_output = 0
    if mceliece == '51220' :
        while anova_matrix[0][0][index_anova_output] < 200:
            index_anova_output += 1
    elif mceliece == '51220' :
        while anova_matrix[0][index_anova_output] < 100:
            index_anova_output += 1
    elif mceliece == '8192128' :
        while (anova_matrix[0][index_anova_output] < 100) or (np.isnan(anova_matrix[0][index_anova_output])):
            index_anova_output += 1
    traces_separation_list.append(index_anova_output)
    print(traces_separation_list[0])

    if mceliece == '51220':
        for index_coefficient in range(1,anova_matrix.shape[1]):
            while (anova_matrix[0][index_coefficient][index_anova_output] < anova_matrix[0][index_coefficient-1][index_anova_output]) or (anova_matrix[0][index_coefficient][index_anova_output] < 100):
                index_anova_output += 1
            traces_separation_list.append(index_anova_output)
            if index_coefficient - 1 : assert traces_separation_list[index_coefficient]-traces_separation_list[index_coefficient-1] == 193
    elif mceliece == '348864' :
        for index_coefficient in range(1,anova_matrix.shape[0]):
            while (anova_matrix[index_coefficient][index_anova_output] < anova_matrix[index_coefficient-1][index_anova_output]) or (anova_matrix[index_coefficient][index_anova_output] < 200):
                index_anova_output += 1
            traces_separation_list.append(index_anova_output)
            print(traces_separation_list[index_coefficient])
            print(traces_separation_list[index_coefficient]-traces_separation_list[index_coefficient-1])
            if index_coefficient - 1 : assert traces_separation_list[index_coefficient]-traces_separation_list[index_coefficient-1] == 251
    elif mceliece == '8192128' :
        for index_coefficient in range(1,anova_matrix.shape[0]):
            while (anova_matrix[index_coefficient][index_anova_output] < anova_matrix[index_coefficient-1][index_anova_output]) or (anova_matrix[index_coefficient][index_anova_output] < 100):
                index_anova_output += 1
            traces_separation_list.append(index_anova_output)
            print(traces_separation_list[index_coefficient])
            print(traces_separation_list[index_coefficient]-traces_separation_list[index_coefficient-1])
            if index_coefficient - 1 : assert traces_separation_list[index_coefficient]-traces_separation_list[index_coefficient-1] == 517

    alpha_traces = np.load(alpha_traces)
    alpha_traces = alpha_traces['traces_matrix']
    print(alpha_traces.shape)
    
    intermediate_values = np.load(intermediate_values)
    intermediate_values = intermediate_values['intermediate_values_matrix']
    
    nb_samples_inv_matrix = traces_separation_list[1] - traces_separation_list[0]
    if mceliece == '8192128' :
        beginning_decapsulation = 0
        nb_decapsulation_traces = 200
        inv_traces_matrix = np.zeros((nb_decapsulation_traces, alpha_traces.shape[1], nb_samples_inv_matrix), dtype="float")
    else : 
        beginning_decapsulation = 0
        nb_decapsulation_traces = alpha_traces.shape[0]
        inv_traces_matrix = np.zeros((alpha_traces.shape[0], alpha_traces.shape[1], nb_samples_inv_matrix), dtype="float")
    print(inv_traces_matrix.shape)

    print(intermediate_values.shape)
    if mceliece == '8192128' : inv_intermediate_values = np.zeros((nb_decapsulation_traces, alpha_traces.shape[1]), dtype=type(intermediate_values[0][0][0]))
    else : inv_intermediate_values = np.zeros((alpha_traces.shape[0], alpha_traces.shape[1]), dtype=type(intermediate_values[0][0][0]))
    print(type(intermediate_values[0][0][0]))
    print(inv_intermediate_values.shape)
    
    nb_samples_loop_matrix = traces_separation_list[2]-traces_separation_list[1]
    if mceliece == '8192128' : loop_traces_matrix = np.zeros((nb_decapsulation_traces, alpha_traces.shape[1], 2*t-1, nb_samples_loop_matrix), dtype="float")
    else : loop_traces_matrix = np.zeros((alpha_traces.shape[0], alpha_traces.shape[1], 2*t-1, nb_samples_loop_matrix), dtype="float")
    print(loop_traces_matrix.shape)

    if mceliece == '8192128' : loop_intermediate_values = np.zeros((nb_decapsulation_traces, alpha_traces.shape[1], 2*t-1), dtype=type(intermediate_values[0][0][0]))
    else : loop_intermediate_values = np.zeros((alpha_traces.shape[0], alpha_traces.shape[1], 2*t-1), dtype=type(intermediate_values[0][0][0]))
    print(loop_intermediate_values.shape)

    for index_trace in tqdm(range(beginning_decapsulation, beginning_decapsulation+nb_decapsulation_traces)):
        for index_alpha in range(alpha_traces.shape[1]):
            if mceliece == '8192128' : 
                inv_traces_matrix[index_trace-beginning_decapsulation][index_alpha][:] = alpha_traces[index_trace][index_alpha][traces_separation_list[0]:traces_separation_list[1]]
                inv_intermediate_values[index_trace-beginning_decapsulation][index_alpha] = intermediate_values[index_trace][index_alpha][0]
            else : 
                inv_traces_matrix[index_trace][index_alpha][:] = alpha_traces[index_trace][index_alpha][traces_separation_list[0]:traces_separation_list[1]]
                inv_intermediate_values[index_trace][index_alpha] = intermediate_values[index_trace][index_alpha][0]
            for index_loop_trace in range(2*t-1):
                if mceliece == '8192128' : 
                    loop_traces_matrix[index_trace-beginning_decapsulation][index_alpha][index_loop_trace][:] = alpha_traces[index_trace][index_alpha][
                                       traces_separation_list[1]+517*index_loop_trace:traces_separation_list[2]+517*index_loop_trace]
                    loop_intermediate_values[index_trace-beginning_decapsulation][index_alpha][index_loop_trace] = intermediate_values[index_trace][index_alpha][index_loop_trace+1]
                else:
                    loop_traces_matrix[index_trace][index_alpha][index_loop_trace][:] = alpha_traces[index_trace][index_alpha][
                                       traces_separation_list[index_loop_trace+1]:traces_separation_list[index_loop_trace+2]]
                    loop_intermediate_values[index_trace][index_alpha][index_loop_trace] = intermediate_values[index_trace][index_alpha][index_loop_trace+1]

    np.savez_compressed(coefficient_intermediate_values, inv_intermediate_values=inv_intermediate_values, loop_intermediate_values=loop_intermediate_values)  
    np.savez_compressed(coefficient_traces, inv_traces_matrix=inv_traces_matrix, loop_traces_matrix=loop_traces_matrix)    

if __name__ == '__main__':
    import argparse
 
    parser = argparse.ArgumentParser()

    parser.add_argument('-l','--leakage_assessment', dest='leakage_assessment', action='store_true', default=False)
    
    parser.add_argument('-d','--decapsulation', action='store_true', dest='decapsulation_slicing', default=False)

    parser.add_argument('-i','--intermediate_values', action='store_true', dest='intermediate_values', default=False)
    
    parser.add_argument('-a','--alpha', action='store_true', dest='alpha_slicing', default=False)
    
    parser.add_argument('--me','--mceliece', action='store', dest='mceliece', required=True,
                        choices=['51220', '348864', '8192128'],help='McEliece security paramaters')

    args = parser.parse_args()
    
    if args.leakage_assessment : 
        parameters_input = 'results/matrices/'+str(args.mceliece)+'/parameters_synd_'+str(args.mceliece)+'.npz'
        list_decapsulation_traces_input = list('results/matrices/'+str(args.mceliece)+'/leakage_assessment/decapsulation_traces_matrix_'+str(args.mceliece)+'_'+str(nb_file)+'.npz' for nb_file in range(10))
        
        alpha_traces_output = 'results/matrices/'+str(args.mceliece)+'/leakage_assessment/alpha_traces_'+str(args.mceliece)+'.npz'
        intermediates_values_output = 'results/matrices/'+str(args.mceliece)+'/leakage_assessment/intermediate_values_'+str(args.mceliece)+'.npz'
    else :
        parameters_input = 'results/matrices/'+str(args.mceliece)+'/parameters_synd_'+str(args.mceliece)+'.npz'
        if args.mceliece == "51220":
            list_decapsulation_traces_input = list('results/matrices/'+str(args.mceliece)+'/attack/decapsulation_traces_matrix_'+str(args.mceliece)+'_attack_'+str(nb_file)+'.npz' for nb_file in range(2))
        elif args.mceliece == "348864":
            list_decapsulation_traces_input = list('results/matrices/'+str(args.mceliece)+'/attack/decapsulation_traces_matrix_'+str(args.mceliece)+'_attack_'+str(nb_file)+'.npz' for nb_file in range(7))
        elif args.mceliece == "8192128":
            beginning_matrix = 4
            nb_matrix_to_load = 2
            list_decapsulation_traces_input = list('results/matrices/'+str(args.mceliece)+'/attack/decapsulation_traces_matrix_'+str(args.mceliece)+'_attack_'+str(nb_file)+'.npz' for nb_file in range(beginning_matrix, beginning_matrix + nb_matrix_to_load))
        
        if args.mceliece == '8192128':
            beginning_decapsulation_matrix = 1
            alpha_traces_output = 'results/matrices/'+str(args.mceliece)+'/attack/alpha_traces_'+str(args.mceliece)+'_'+str(beginning_decapsulation_matrix)+'.npz'
            intermediates_values_output = 'results/matrices/'+str(args.mceliece)+'/attack/intermediate_values_'+str(args.mceliece)+'_'+str(beginning_decapsulation_matrix)+'.npz'
            nb_coefficient_matrix = 2
            coefficient_traces = 'results/matrices/'+str(args.mceliece)+'/attack/coefficient_traces_'+str(args.mceliece)+'_'+str(nb_coefficient_matrix)+'.npz'
            coefficient_intermediates_values = 'results/matrices/'+str(args.mceliece)+'/attack/coefficient_intermediates_values_'+str(args.mceliece)+'_'+str(nb_coefficient_matrix)+'.npz'  
        else:
            alpha_traces_output = 'results/matrices/'+str(args.mceliece)+'/attack/alpha_traces_'+str(args.mceliece)+'.npz'
            intermediates_values_output = 'results/matrices/'+str(args.mceliece)+'/attack/intermediate_values_'+str(args.mceliece)+'.npz'
            coefficient_traces = 'results/matrices/'+str(args.mceliece)+'/attack/coefficient_traces_'+str(args.mceliece)+'.npz'
            coefficient_intermediates_values = 'results/matrices/'+str(args.mceliece)+'/attack/coefficient_intermediates_values_'+str(args.mceliece)+'.npz'

    anova_test_output = 'results/statistical_tests/'+str(args.mceliece)+'/anova_test_'+str(args.mceliece)+'.npz'

    print(list_decapsulation_traces_input)
    print(coefficient_traces)
    if args.decapsulation_slicing : decapsulation_traces_slicing(alpha_traces_output, list_decapsulation_traces_input, args.mceliece)
    if args.intermediate_values : intermediate_values_computation(intermediates_values_output, alpha_traces_output, parameters_input, args.mceliece)
    if args.alpha_slicing : alpha_traces_slicing(coefficient_traces, coefficient_intermediates_values, anova_test_output, alpha_traces_output, intermediates_values_output, args.mceliece)