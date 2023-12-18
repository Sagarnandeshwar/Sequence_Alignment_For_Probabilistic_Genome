import numpy as np
import random
import time

file1 = "./dataset/chr22.maf.ancestors.42000000.complete.boreo.fa.txt"
file2 = "./dataset/chr22.maf.ancestors.42000000.complete.boreo.conf.txt"


def readData(filename):
    with open(filename) as f:
        data = f.readline()
    return data


def processData(ch22, prob_float):
    dataset = np.empty((2,len(ch22)),dtype = object)

    for ind in range(len(ch22)):
        dataset[0][ind] = ch22[ind]
        dataset[1][ind] = float(prob_float[ind])

    return dataset


def getData(f1, f2):
    ch22 = readData(f1)
    prob = readData(f2)
    prob = prob.split()
    prob_float = []
    for p in prob:
        prob_float.append(p)

    dataset = processData(ch22, prob_float)
    return ch22, prob_float, dataset


def list_to_str(lis):
    str = ""
    for element in lis:
        str = str + element
    return str


def generatQuery(n, query_length, genome, indel, subs):
    queries = []
    starts = np.random.randint(0, high=len(genome)-1-query_length, size=n)
    for s in starts:
        queries.append(genome[s:s+query_length])

    queries_insert = []
    for query in queries:
        new_query = []
        ind = 0
        while ind < len(query):
            insert = np.random.binomial(1, indel)
            if insert:
                nuc = random.choices(['A', 'C', 'G', 'T'])[0]
                new_query.append(nuc)
                continue
            new_query.append(query[ind])
            ind = ind + 1
        queries_insert.append(list_to_str(new_query))


    queries_delete = []
    for query in queries_insert:
        new_query = []
        ind = 0
        while ind < len(query):
            delete = np.random.binomial(1, indel)
            if delete:
                ind = ind + 1
                continue
            new_query.append(query[ind])
            ind = ind + 1
        queries_delete.append(list_to_str(new_query))

    queries_sub = []
    for query in queries_delete:
        new_query = []
        ind = 0
        while ind < len(query):
            substitute = np.random.binomial(1,subs)
            if substitute:
                nuc = random.choices(['A','C','G','T'])[0]
                new_query.append(nuc)
                ind = ind + 1
                continue
            new_query.append(query[ind])
            ind = ind + 1
        queries_sub.append(list_to_str(new_query))

    """print(queries)
    print(queries_insert)
    print(queries_delete)
    print(queries_sub)"""

    return queries_sub, starts


def findPerfectMatch(index_start_end, query, genome, w_size):
    word = query[index_start_end[0]:index_start_end[1]]
    matches = []
    for ind in range(len(genome)-w_size):
        if genome[ind:ind+w_size-1] == word:
            matches.append([ind, ind+w_size-1])
    return matches


def cal_score(qs, qe, gs, ge, query, genome, process_data):
    length = qe - qs
    ind_q = qs
    ind_g = gs
    score = 0
    while ind_q < length:
        if genome[ind_g] == query[ind_q]:
            score = score + process_data[1][ind_g]
        else:
            score = score - process_data[1][ind_g]
        ind_q = ind_q + 1
        ind_g = ind_g + 1
    return score


def blast(query, genome, process_data, w_size = 11):
    w_mer_index = []
    for ind in range(len(query)-w_size):
        w_mer_index.append((ind,ind+w_size-1))

    seeds = {}
    for index_start_end in w_mer_index:
        matches = findPerfectMatch(index_start_end, query, genome, w_size)
        if len(matches) != 0:
            seeds[index_start_end] = matches

    print(len(seeds))

    seed_DB = {}
    for query_index in seeds.keys():
        list_genome_index = seeds[query_index]
        for genome_index in list_genome_index:
            score = cal_score(query_index[0], query_index[1], genome_index[0], genome_index[1], query, genome, process_data)
            seed_DB[(query_index[0], query_index[1], genome_index[0], genome_index[1])] = score

    ungap_extension = {}
    for index in seed_DB.keys():
        score_index = unGapExtension(index, query, genome, process_data, seed_DB[index])
        ungap_extension[index] = score_index

    ungap_extension_DB = {}
    for index in ungap_extension.keys():
        [s, qs, qe, gs, ge] = ungap_extension[index]
        if s > threshold2:
            ungap_extension_DB[(qs, qe, gs, ge)] = s

    gap_extension_DB = {}
    for index in ungap_extension_DB.keys():
        (s, al_q, al_g, qs, qe, gs, ge) = gapExtension(index, query, genome, process_data, ungap_extension_DB[index])
        gap_extension_DB[(qs, qe, gs, ge)] = {"score": s, "query": al_q, "genome": al_g}

    return gap_extension_DB


def unGapExtension(index, query, genome, process_data, prevScore):
    Q_END = len(query)
    G_END = len(genome)
    query_start = index[0]
    query_end = index[1]
    genome_start = index[2]
    genome_end = index[3]
    leftExtension = True
    rightExtention = True

    left_score = 0
    right_score = 0
    max_left_score = 0
    max_right_score = 0

    high_left_query = query_start
    high_left_genome = genome_start
    high_right_query = query_end
    high_right_genome = genome_end

    while leftExtension or rightExtention:
        if rightExtention:
            if (query_end + 1 > Q_END - 1) or (genome_end + 1 > G_END - 1) :
                rightExtention = False
                continue
            query_end = query_end + 1
            genome_end = genome_end + 1
            nuc_genome = genome[genome_end]
            right_pos_query = query[query_end]
            if nuc_genome == right_pos_query:
                right_score = right_score + process_data[1][genome_end] - (1-process_data[1][genome_end])
            else:
                right_score = right_score - process_data[1][genome_end] + (process_data[1][genome_end]*1/3) - (process_data[1][genome_end]*2/3)

            if right_score >= max_right_score:
                max_right_score = right_score
                high_right_query = query_end
                high_right_genome = genome_end

            if right_score + threshold1 < max_right_score:
                rightExtention = False

        if leftExtension:
            if query_start - 1 < 0 or genome_start - 1 < 0:
                leftExtension = False
                continue
            query_start = query_start - 1
            genome_start = genome_start - 1
            nuc_genome = genome[query_start]
            left_pos_query = query[query_start]
            if nuc_genome == left_pos_query:
                left_score = left_score + process_data[1][genome_start] - (1 - process_data[1][genome_start])
            else:
                left_score = left_score - process_data[1][genome_start] + (
                            (1 - process_data[1][genome_start]) * 1 / 3) - ((1 - process_data[1][genome_start]) * 2 / 3)

            if left_score >= max_left_score:
                max_left_score = left_score
                high_left_query = query_start
                high_left_genome = genome_start

            if left_score + threshold1 < max_left_score:
                leftExtension = False

    score = prevScore + max_right_score + max_left_score
    return [score, high_left_query, high_right_query, high_left_genome, high_right_genome]


def gapExtension(index, query, genome, process_data, prev_score):
    """
    0 -> diagnol
    1 -> up
    2 -> left
    """

    Q_END = len(query) - 1
    G_END = len(genome) - 1

    query_left = ""
    genome_left = ""
    query_right = ""
    genome_right = ""
    query_start = index[0]
    query_end = index[1]
    genome_start = index[2]
    genome_end = index[3]
    left_max_score = 0
    right_max_score = 0

    high_left_query = query_start
    high_left_genome = genome_start

    high_right_query = query_end
    high_right_genome = genome_end

    if query_end < len(query) - 1:
        dp_matrix = np.zeros((len(genome) - genome_end - 1, len(query) - query_end), dtype=float)
        path = np.zeros((len(genome) - genome_end - 1, len(query) - query_end), dtype=int)
        row_length, col_length = dp_matrix.shape
        for row in range(1, row_length):
            dp_matrix[row][0] = dp_matrix[row - 1][0] - 1
            path[row][0] = 1
        for col in range(1, col_length):
            dp_matrix[0][col] = dp_matrix[0][col - 1] - 1
            path[0][col] = 2
        for row in range(1, row_length):
            for col in range(1, col_length):
                diagonal = dp_matrix[row - 1][col - 1]
                if genome[genome_end + row] == query[query_end + col]:
                    diagonal = diagonal + process_data[1][genome_end + row] - (1 - process_data[1][genome_end + row])
                else:
                    diagonal = diagonal - process_data[1][genome_end + row] + (1 - process_data[1][genome_end + row])
                up = dp_matrix[row - 1][col] - 1
                left = dp_matrix[row][col - 1] - 1
                if diagonal > up:
                    if diagonal > left:
                        dp_matrix[row][col] = diagonal
                        path[row][col] = 0
                    else:
                        dp_matrix[row][col] = left
                        path[row][col] = 2
                else:
                    if up > left:
                        dp_matrix[row][col] = up
                        path[row][col] = 1
                    else:
                        dp_matrix[row][col] = left
                        path[row][col] = 2

        last = [row[-1] for row in dp_matrix]
        right_max_score = max(last)
        col = len(dp_matrix[0]) - 1
        row = last.index(right_max_score)
        high_right_query = high_right_query + col
        high_right_genome = high_right_genome + row
        while row > 0 and col > 0:
            if path[row][col] == 0:
                query_right = query[query_end + col] + query_right
                genome_right = genome[genome_end + row] + genome_right
                row = row - 1
                col = col - 1
            elif path[row][col] == 1:
                query_right = "-" + query_right
                genome_right = genome[genome_end + row] + genome_right
                row = row - 1
            else:
                query_right = query[query_end + col] + query_right
                genome_right = "-" + genome_right
                col = col - 1

    if query_start > 0:
        dp_matrix = np.zeros((genome_start + 1, query_start + 1), dtype=float)
        path = np.zeros((genome_start + 1, query_start + 1), dtype=int)
        row_length, col_length = dp_matrix.shape
        for row in range(1, row_length):
            dp_matrix[row][0] = dp_matrix[row - 1][0] - 1

            path[row][0] = 1
        for col in range(1, col_length):
            dp_matrix[0][col] = dp_matrix[0][col - 1] - 1
            path[0][col] = 2
        for row in range(1, row_length):
            for col in range(1, col_length):
                diagonal = dp_matrix[row - 1][col - 1]
                if genome[genome_start - row] == query[query_start - col]:
                    diagonal = diagonal + process_data[1][genome_start - row] - (1 - process_data[1][genome_start - row])
                else:
                    diagonal = diagonal - process_data[1][genome_start - row] + (1 - process_data[1][genome_start - row])
                up = dp_matrix[row - 1][col] - 1
                left = dp_matrix[row][col - 1] - 1
                if diagonal > up:
                    if diagonal > left:
                        dp_matrix[row][col] = diagonal
                        path[row][col] = 0
                    else:
                        dp_matrix[row][col] = left
                        path[row][col] = 2
                else:
                    if up > left:
                        dp_matrix[row][col] = up
                        path[row][col] = 1
                    else:
                        dp_matrix[row][col] = left
                        path[row][col] = 2

        last_col = [row[-1] for row in dp_matrix]
        left_max_score = max(last_col)
        col = len(dp_matrix[0]) - 1
        row = last_col.index(left_max_score)
        high_left_query = high_left_query - col
        high_left_genome = high_left_genome - row
        while row > 0 and col > 0:
            if path[row][col] == 0:
                query_left = query_left + query[query_start - col]
                genome_left = genome_left + genome[genome_start - row]
                row = row - 1
                col = col - 1
            elif path[row][col] == 1:
                query_left = query_left + "-"
                genome_left = genome_left + genome[genome_start - row]
                row = row - 1
            else:
                query_left = query_left + query[query_start - col]
                genome_left = genome_left + "-"
                col = col - 1

    align_query = query_left + query[query_start:query_end] + query_right
    align_genome = genome_left + genome[genome_start:genome_end] + genome_right
    score = prev_score + left_max_score + right_max_score

    """print("left nw")
    print(query_left)
    print(genome_left)
    print("rigth nw")
    print(query_right)
    print(genome_right)
    print("original query")
    print(query[query_start:query_end])
    print(genome[genome_start:genome_end])
    print("align solution")
    print(align_query)
    print(align_genome)
    print("Final Score")
    print(score)"""
    return (score, align_query, align_genome, high_left_query, high_right_query, high_left_genome, high_right_genome)

# Parameters
threshold1 = 3
threshold2 = 5
indels = 0.10
substitution = 0.10

if __name__ == "__main__":
    genome, prob, process_data = getData(file1, file2)
    queries, start = generatQuery(1, 100, genome, indels, substitution)
    print(start)
    start_time = time.time()
    for a_query in queries:
        print("START")
        print(a_query)
        infoTable = blast(a_query, genome, process_data, 11)
        for index in infoTable.keys():
            print(index)
            print(infoTable[index]["score"])
            print(infoTable[index]["query"])
            print(infoTable[index]["genome"])
        print("END")
        print("--- %s seconds ---" % (time.time() - start_time))