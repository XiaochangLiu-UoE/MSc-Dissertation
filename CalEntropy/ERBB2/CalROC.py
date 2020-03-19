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


def count_ap(mols):
    count = 0
    for mol in mols:
        props = mol[1]
        if props == 1:
            count = count + 1
    return count


def read_ranking(output):
    ranking = {}
    ligand = []
    with open(output, 'r') as o:
        for line in o:
            line = line.rstrip()
            if line == '> <Name>':
                name = next(o).rstrip()
            if line == '>  <SCORE_INFO>':
                score = next(o).rstrip().split(' ')
                ligand.append(float(score[0]))
            if line == '$$$$':
                if name not in ranking:
                    if 'ZINC' in name:
                        ligand.append(0)
                    elif 'BDB' in name:
                        ligand.append(1)
                    ranking[name] = ligand
                ligand = []
    shortlist = list(ranking.values())
    num_ap = count_ap(mols=shortlist)
    n1 = 40 - num_ap
    n2 = 1240 - len(shortlist) - n1
    if n1 != 0:
        for i in range(n1):
            shortlist.append([1, 1])
    if n2 != 0:
        for i in range(n2):
            shortlist.append([1, 0])
    return shortlist


output = 'output.sdf'
tpr, fpr, auc = calculate_auc_roc(ranking=read_ranking(output=output))
print(auc)