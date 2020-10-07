class GetData:
    @staticmethod
    def get_reference_data(short_path, doc_name):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open(short_path + doc_name, 'r') as ref_data:
            for line in ref_data:
                line = line.strip().split('\t')

                # other loci - there is int meaning of lr
                if len(line) == 4:
                    locus = line[0]
                    lr = line[3]

                    #  case of gender specific loci
                    if lr == '-':
                        ref_dict[locus] = lr
                        continue
                    # case of int lr meanings of loci
                    else:
                        lr = float(line[3]) * 100 / 100
                        ref_dict[locus] = lr

                # case for getting cpi and p meanings
                elif len(line) == 1 and line[0] != '':
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    elif line[0] == 'P':
                        p = float(line[1])
        return ref_dict, cpi, p