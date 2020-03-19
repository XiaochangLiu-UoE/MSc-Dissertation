import csv


def calculate_auc_roc(ranking):
    y = []
    for rank in ranking:
        y.append(rank[1])
    y = tuple(y)
    tpr = []
    fpr = []
    sfp = 0
    p = 0
    ap = y.count(1)
    an = len(y) - ap
    tp = 0.0
    fp = 0.0
    for j in range(len(ranking)):
        if j != len(ranking) - 1 and ranking[j][0] == ranking[j + 1][0]:
            if ranking[j][1] == 0:
                sfp = sfp + 1
            true_positive_rate = tp / ap
            false_positive_rate = fp / an
            tpr.append(true_positive_rate)
            fpr.append(false_positive_rate)
        else:
            p = j + 1
            if ranking[j][1] == 0:
                fp = fp + 1
                sfp = sfp + 1
            fp = sfp
            tp = p - fp
            true_positive_rate = tp / ap
            false_positive_rate = fp / an
            tpr.append(true_positive_rate)
            fpr.append(false_positive_rate)
    tpr = tuple(tpr)
    fpr = tuple(fpr)
    area = [[0.0, 0.0]]
    for i in range(len(fpr)):
        area.append(tuple([fpr[i], tpr[i]]))
    area = tuple(area)
    auc = 0.0
    start_point = 0.0
    end_point = 0.0
    height = 0.0
    for key in area:
        if start_point != key[0] and height != key[1]:
            end_point = key[0]
            temp = ((end_point - start_point) * (height + key[1])) / 2
            auc = auc + temp
            start_point = key[0]
            height = key[1]
        elif start_point != key[0] or height != key[1]:
            end_point = key[0]
            temp = (end_point - start_point) * height
            auc = auc + temp
            start_point = key[0]
            height = key[1]
    return tpr, fpr, auc


def tell_active(name):
    if 'BDB' in name:
        return 1
    if 'ZINC' in name:
        return 0


def read_csv(csvfile, heading=True):
    if heading:
        with open(csvfile, 'r') as csv:
            lines = tuple(tuple(line.strip().split(',')) for line in csv.readlines()[1:])
        return lines
    else:
        with open(csvfile, 'r') as csv:
            lines = tuple(tuple(line.strip().split(',')) for line in csv.readlines())
    return lines


def return_ranking(output):
    output = list(sorted(iter(output), key=lambda m: m[1], reverse=False))
    ranking = {}
    lds = []
    ap_count = 0
    for mol in output:
        lds.append(mol[1])
        if mol[0] in ranking:
            lds = []
        else:
            active = tell_active(mol[0])
            lds.append(active)
            if active == 1:
                ap_count = ap_count + 1
            ranking[mol[0]] = lds
            lds = []
    true_rank = list(ranking.values())
    left_actives = 40 - ap_count
    left_decoys = 1240 - left_actives - len(true_rank)
    if left_actives != 0:
        for i in range(left_actives):
            true_rank.append((9999, 1))
    if left_decoys != 0:
        for i in range(left_decoys):
            true_rank.append((9999, 0))
    return tuple(true_rank)


def return_output(mols):
    output = []
    for mol in mols:
        output.append((mol[0], float(mol[1]) + float(mol[5])))
    return tuple(output)


csvfile = 'Entropy.csv'
mols = read_csv(csvfile=csvfile)
output = return_output(mols=mols)
ranking = return_ranking(output=output)
tpr, fpr, auc = calculate_auc_roc(ranking=ranking)
print(auc)