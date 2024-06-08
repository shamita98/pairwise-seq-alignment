# function to perform backward sequence alignment tracing
def seq_tracing(seq_i, seq_j, route_matrix, new_seq_i, new_seq_j, current_loc_i, current_loc_j):
    current_i = current_loc_i
    current_j = current_loc_j
    previous_i = route_matrix[current_i][current_j][0]
    previous_j = route_matrix[current_i][current_j][1]
    if previous_i + 1 == current_i and previous_j == current_j:  # down-path
        new_seq_i += '-'
        new_seq_j += seq_j[previous_i]
    elif previous_i + 1 == current_i and previous_j + 1 == current_j:  # diagonal-path
        new_seq_i += seq_i[previous_j]
        new_seq_j += seq_j[previous_i]
    elif previous_i == current_i and previous_j + 1 == current_j:  # right-path
        new_seq_i += seq_i[previous_j]
        new_seq_j += '-'
    current_i = previous_i
    current_j = previous_j
    return current_i, current_j, new_seq_i, new_seq_j


# function to reverse the aligned sequence after tracing
def reverse(string):
    length = len(string)
    i = length - 1
    new_string = ''
    while i > 0 or i == 0:
        new_string += string[i]
        i -= 1
    return new_string


# enter sequences to be aligned, penalty, match and mismatch scores
seq_i = input('Enter sequence 1: ').upper()
seq_j = input('Enter sequence 2: ').upper()
gap_penalty = int(input('Enter gap penalty score: '))
match = int(input('Enter match score: '))
mismatch = int(input('Enter mismatch score: '))


# create a score matrix with row 0 and col 0 filled with penalty scores
num_i = len(seq_j) + 1  # row
num_j = len(seq_i) + 1  # col
score_matrix = []
route_matrix = []
for i in range(num_i):
    score_matrix.append(num_j * [0])  # num_j*[0] = columns
    row_in_route_matrix = []
    if i == 0:
        row_in_route_matrix.append([-1, -1])
        for j in range(1, num_j):
            score_matrix[i][j] = gap_penalty + gap_penalty * (j - 1)
            row_in_route_matrix.append([i, j - 1])
    elif i > 0:
        score_matrix[i][0] = gap_penalty + gap_penalty * (i - 1)
        row_in_route_matrix.append([i - 1, 0])
    route_matrix.append(row_in_route_matrix)


# fill each box with score considering all the 3 possible paths (diagonal, down and right)
for i in range(1, num_i):
    for j in range(1, num_j):
        if seq_j[i - 1] == seq_i[j - 1]:
            diagonal_score = score_matrix[i - 1][j - 1] + match
        else:
            diagonal_score = score_matrix[i - 1][j - 1] + mismatch
            # diagonal_score = score_matrix[i - 1][j - 1] + sub_matrix[seq_j[i - 1]][seq_i[j - 1]]
        right_score = score_matrix[i][j - 1] + gap_penalty
        down_score = score_matrix[i - 1][j] + gap_penalty
        previous_score_loc = [[i - 1, j - 1], [i, j - 1], [i - 1, j]]
        best_score = max(diagonal_score, right_score, down_score)
        score_matrix[i][j] = best_score

        max_total_score = best_score + score_matrix[i - 1][j - 1]  # filling in route matrix
        best_previous_score_loc = [i - 1, j - 1]
        for loc in previous_score_loc:
            total_score = best_score + score_matrix[loc[0]][loc[1]]
            if total_score > max_total_score:
                max_total_score = total_score
                best_previous_score_loc = loc
        route_matrix[i].append(best_previous_score_loc)


# building alignment by tracing from the lowest right corner using seq_tracing()
current_i = num_i - 1
current_j = num_j - 1
new_seq_i = ''
new_seq_j = ''

while current_i > 0 or current_j > 0:
    if current_i == 0 and current_j == 0:
        break
    else:
        current_i, current_j, new_seq_i, new_seq_j = seq_tracing(seq_i, seq_j, route_matrix,
                                                               new_seq_i, new_seq_j, current_i, current_j)

new_seq_i = reverse(new_seq_i)
new_seq_j = reverse(new_seq_j)
annotation_line = ''
for i in range(len(new_seq_j)):
    if new_seq_i[i] == new_seq_j[i]:
        annotation_line += '|'
    elif new_seq_i[i] == '-' or new_seq_j[i] == '-':
        annotation_line += ' '
    else:
        annotation_line += '*'

print(f'{new_seq_i}\n{annotation_line}\n{new_seq_j}')








