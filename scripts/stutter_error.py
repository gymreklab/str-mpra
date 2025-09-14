import numpy as np
import pandas as pd
from collections import defaultdict
import csv
import argparse

class StutterAnalyzer:
    def __init__(self, file_path, out_path=None):
        self.df = pd.read_csv(file_path, sep='\t', names=["Reads number", "STR_name", "BC"])
        self.result_dict = {}
        self.dup_BC = {}
        self.valid_dup_BC = {}
        self.invalid_dup_BC = {}
        self.uniq_BC = {}
        self.readNum_match = 0
        self.merged_BC = {}
        self.out_path = out_path

        self.analyze_data()

    def analyze_data(self):
        for _, row in self.df.iterrows():
            BC = row["BC"]
            name = row["STR_name"]
            reads_num = row["Reads number"]

            if BC not in self.result_dict:
                self.result_dict[BC] = []

            self.result_dict[BC].append(f"{name}_{reads_num}")
       
        for key, value in self.result_dict.items():
            if len(value) >= 2:
                self.dup_BC[key] = value
            else:
                self.readNum_match += int(value[0].split("_")[-1])
                self.uniq_BC[key] = value[0]

        for key, value in self.dup_BC.items():
            STR = None
            same_group = True

            for v in value:
                if STR is None:
                    STR = v.split("_")[:-3]
                elif STR != v.split("_")[:-3]:
                    same_group = False
                    break

            if same_group:
                rpts = []
                reads_num = []

                for v in value:
                    name = "_".join(v.split("_")[:-3])
                    motif = v.split("_")[-3]
                    rpts.append(v.split("_")[-2])
                    reads_num.append(v.split("_")[-1])

                data = {
                    "STR_name": name,
                    "motif": motif,
                    "motif_length": len(motif),
                    "rpt": rpts,
                    "reads_num": reads_num
                }

                self.valid_dup_BC[key] = data
            else:
                self.invalid_dup_BC[key] = value

    def bestGuess(self, u, d, iterations=100, threshold=0.75):
        for key, value in self.valid_dup_BC.items():
            self.valid_dup_BC[key] = StutterAnalyzer.possible_rpt(value, u, d,)

        result = self.updateParams()
        new_params = result[:2] 

        for i in range(iterations):  
            for key, value in self.valid_dup_BC.items():
                self.valid_dup_BC[key] = StutterAnalyzer.possible_rpt(value, *new_params,)
            
            current_params = new_params
            result = self.updateParams()
            new_params = result[:2] 
            
            delta = tuple(x - y for x, y in zip(current_params, new_params))
            
            if all(abs(x) < 0.001 for x in delta[:2]): 
                break
        
        for k, v in self.valid_dup_BC.items():
            if v["best_guess"]["QS"] >= threshold:
                joined_str = '_'.join([
                    v["STR_name"],
                    v["motif"],
                    v["best_guess"]["best_rpt"],
                    str(sum(list(map(int, v["reads_num"])))),
                ])
                self.merged_BC[k] = joined_str
        
        self.merged_BC.update(self.uniq_BC)
       
        if self.out_path is None:
            print("Output path is not specified, no output is saved")
        else:
            with open(self.out_path, 'w', newline='') as tsvfile:
                writer = csv.writer(tsvfile, delimiter='\t')
                for k,v in self.merged_BC.items():
                    v = v.rsplit('_', 1)
                    writer.writerow([v[1], k, v[0]])

        return delta, new_params, i, result[2]


    def updateParams(self):
        delta_r_sum = defaultdict(int)

        for k, v in self.valid_dup_BC.items():

            delta_rs = v['best_guess']['best_delta_r']
            reads_num = v['reads_num']

            for r, n in zip(delta_rs, reads_num):
                delta_r_sum[r] += int(n)

        ref_case = delete = add = 0
        sum_x_u = 0
        sum_x_d = 0
        for k, v in delta_r_sum.items():
            if int(k) == 0:
                ref_case = v
            elif int(k) < 0:
                delete += v
                sum_x_d += -k*v
            else:
                add += v
                sum_x_u += k*v
        s = ref_case + delete + add + self.readNum_match
 
        return add/s, delete/s, delta_r_sum

    @staticmethod 
    def convert_num(input_array):
        output = []
        for i in input_array:
            if i == "ref":
                output.append(0)
            elif i[0] == 'm':
                output.append(int("-" + i[1:]))
            elif i[0] == 'p':
                output.append(int("+" + i[1:]))
        return output

    @staticmethod
    def stutter_error_nonGeo_log_prob(delta_r, u, d):
        if delta_r == 0:
            return np.log(1 - u - d)
        elif delta_r > 0:
            return np.log(u)
        else:
            return np.log(d)

    @staticmethod
    def relative_delta_r(input_array, index):
        return np.array(input_array, dtype = int) - int(input_array[index])

    @staticmethod
    def loglikelihood(a, u, d,):
        rpt_num = StutterAnalyzer.convert_num(a["rpt"])

        LogLikelihood = []

        for i in range(len(rpt_num)):
            r_array = StutterAnalyzer.relative_delta_r(rpt_num, i)
            ll = 0
            for j in range(len(r_array)):
                se_prob = StutterAnalyzer.stutter_error_nonGeo_log_prob(r_array[j], u, d)
                occurrence = int(a["reads_num"][j])
                ll += se_prob * occurrence

            LogLikelihood.append(ll)
        return LogLikelihood
    
    @staticmethod
    def possible_rpt(a, u, d,):
        loglikelihoods = StutterAnalyzer.loglikelihood(a, u, d,)  
        max_index = loglikelihoods.index(max(loglikelihoods))
        quality_score = np.exp(max(loglikelihoods))/sum(np.exp(loglikelihoods))
        
        rpt_num = StutterAnalyzer.convert_num(a["rpt"])
        best_delta_r = StutterAnalyzer.relative_delta_r(rpt_num, max_index)

        output_dict = a.copy()
        output_dict["likelihoods"] = np.exp(loglikelihoods)
        output_dict["best_guess"] = {
            "QS": quality_score, 
            "best_rpt": a["rpt"][max_index], 
            "best_delta_r": list(best_delta_r),
            }

        return output_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Stutter Analyzer Script')
    parser.add_argument('-i', '--input_file_path', type=str, required=True, help='Input association table before stutter error correction')
    parser.add_argument('-o', '--output_file_path', type=str, required=True, help='Output association table after stutter error correction')
    parser.add_argument('-u', '--initial_u', type=float, default=0.05, help='Initial setting of u')
    parser.add_argument('-d', '--initial_d', type=float, default=0.05, help='Initial setting of d')
    parser.add_argument('-n', '--iterations', type=int, default=100, help='The number of max iterations, default is 100')
    parser.add_argument('-t', '--threshold', type=float, default=0.75, help='Threshold of the QS to select valid duplicated BCs, default is 0.75')

    args = parser.parse_args()

    str_analyzer = StutterAnalyzer(args.input_file_path, args.output_file_path)
    print(f"Input file is {args.input_file_path}")
    result = str_analyzer.bestGuess(args.initial_u, args.initial_d, iterations=args.iterations, threshold=args.threshold)
    print(f"The stutter error corrected association table is saved in {args.output_file_path}")

    