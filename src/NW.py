from numpy import full
import pandas as pd

#build an array of zeroes
seq1 = "ATGGT"
seq2 = "ATGGAT"
# seq1 = "GCTGATTC"
# seq2 = "GATCTGATTA"
n_rows = len("-"+seq1)
n_columns = len("-"+seq2)

scoring_array = full([n_rows,n_columns],0)
print("Scoring array:\n",scoring_array)

traceback_array = full([n_rows,n_columns],"-")
print("Traceback array:\n",traceback_array)

def pretty_table_from_array(data_array, row_labels,col_labels):
    """Show an HTML table from a 2d numpy array"""
    df = pd.DataFrame(data_array,index=row_labels,columns=col_labels)
    print(df)
    #table_html = df.to_html()
    #return HTML(table_html)
row_labels = [label for label in "-"+seq1]
column_labels = [label for label in "-"+seq2]

print("Scoring array:")
pretty_table_from_array(scoring_array,row_labels,column_labels)
print("Traceback array:")
pretty_table_from_array(traceback_array,row_labels,column_labels)

count = 0
for row_index in range(n_rows):
    for col_index in range(n_columns):
        scoring_array[row_index, col_index] = count
        count += 1

pretty_table_from_array(scoring_array, row_labels, column_labels)

up_arrow = "\u2191"
right_arrow = "\u2192"
down_arrow = "\u2193"
left_arrow = "\u2190"
down_right_arrow = "\u2198"
up_left_arrow = "\u2196"

print("Up arrow",up_arrow)
print("Left arrow",left_arrow)
print("Up Left arrow",up_left_arrow)

# build an array of zeroes
n_rows = len(seq1) + 1  # need an extra row up top
n_columns = len(seq2) + 1  # need an extra column on the left
row_labels = [label for label in "-" + seq1]
column_labels = [label for label in "-" + seq2]

scoring_array = full([n_rows, n_columns], 0)
traceback_array = full([n_rows, n_columns], "-")

# Define Unicode arrows we'll use in the traceback array
up_arrow = "\u2191"
right_arrow = "\u2192"
down_arrow = "\u2193"
left_arrow = "\u2190"
down_right_arrow = "\u2198"
up_left_arrow = "\u2196"

arrow = "-"
gap_penalty = -1
match_bonus = 2
mismatch_penalty = -1
# iterate over columns first because we want to do
# all the columns for row 1 before row 2
for row in range(n_rows):
    for col in range(n_columns):
        if row == 0 and col == 0:
            # We're in the upper right corner
            score = 0
            arrow = "-"
        elif row == 0:
            # We're on the first row
            # but NOT in the corner

            # Look up the score of the previous cell (to the left) in the score array\
            previous_score = scoring_array[row, col - 1]
            # add the gap penalty to it's score
            score = previous_score + gap_penalty
            arrow = left_arrow
        elif col == 0:
            # We're on the first column but not in the first row
            previous_score = scoring_array[row - 1, col]
            score = previous_score + gap_penalty
            arrow = up_arrow
        else:
            # We're in a 'middle' cell of the alignment

            # Calculate the scores for coming from above,
            # from the left, (representing an insertion into seq1)
            cell_to_the_left = scoring_array[row, col - 1]
            from_left_score = cell_to_the_left + gap_penalty

            # or from above (representing an insertion into seq2)
            above_cell = scoring_array[row - 1, col]
            from_above_score = above_cell + gap_penalty

            # diagonal cell, representing a substitution (e.g. A --> T)
            diagonal_left_cell = scoring_array[row - 1, col - 1]

            # NOTE: since the table has an extra row and column (the blank ones),
            # when indexing back to the sequence we want row -1 and col - 1.
            # since row 1 represents character 0 of the sequence.
            if seq1[row - 1] == seq2[col - 1]:
                diagonal_left_cell_score = diagonal_left_cell + match_bonus
            else:
                diagonal_left_cell_score = diagonal_left_cell + mismatch_penalty

            score = max([from_left_score, from_above_score, diagonal_left_cell_score])
            # take the max

            # make note of which cell was the max in the traceback array
            # using Unicode arrows
            if score == from_left_score:
                arrow = left_arrow
            elif score == from_above_score:
                arrow = up_arrow
            elif score == diagonal_left_cell_score:
                arrow = up_left_arrow

        traceback_array[row, col] = arrow
        scoring_array[row, col] = score

pretty_table_from_array(scoring_array, row_labels, column_labels)
pretty_table_from_array(traceback_array, row_labels, column_labels)


def traceback_alignment(traceback_array, seq1, seq2, up_arrow="\u2191", \
                        left_arrow="\u2190", up_left_arrow="\u2196", stop="-"):
    """Align seq1 and seq2 using the traceback matrix and return as two strings

    traceback_array -- a numpy array with arrow characters indicating the direction from
    which the best path to a given alignment position originated

    seq1 - a sequence represented as a string
    seq2 - a sequence represented as a string
    up_arrow - the unicode used for the up arrows (there are several arrow symbols in Unicode)
    left_arrow - the unicode used for the left arrows
    up_left_arrow - the unicode used for the diagonal arrows
    stop - the symbol used in the upper left to indicate the end of the alignment
    """

    n_rows = len(seq1) + 1  # need an extra row up top
    n_columns = len(seq2) + 1  # need an extra row up top

    row = len(seq1)
    col = len(seq2)
    arrow = traceback_array[row, col]
    aligned_seq1 = ""
    aligned_seq2 = ""
    alignment_indicator = ""
    while arrow != "-":
        print("Currently on row:", row)
        print("Currently on col:", col)
        arrow = traceback_array[row, col]
        print("Arrow:", arrow)

        if arrow == up_arrow:
            print("insert indel into top sequence")
            # We want to add the new indel onto the left
            # side of the growing aligned sequence
            aligned_seq2 = "-" + aligned_seq2
            aligned_seq1 = seq1[row - 1] + aligned_seq1
            alignment_indicator = " " + alignment_indicator
            row -= 1

        elif arrow == up_left_arrow:
            print("match or mismatch")
            # Note that we look up the row-1 and col-1 indexes
            # because there is an extra "-" character at the
            # start of each sequence
            seq1_character = seq1[row - 1]
            seq2_character = seq2[col - 1]
            aligned_seq1 = seq1[row - 1] + aligned_seq1
            aligned_seq2 = seq2[col - 1] + aligned_seq2
            if seq1_character == seq2_character:
                alignment_indicator = "|" + alignment_indicator
            else:
                alignment_indicator = " " + alignment_indicator
            row -= 1
            col -= 1

        elif arrow == left_arrow:
            print("Insert indel into left sequence")
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[col - 1] + aligned_seq2
            alignment_indicator = " " + alignment_indicator
            col -= 1

        elif arrow == stop:
            break
        else:
            raise ValueError(
                f"Traceback array entry at {row},{col}: {arrow} is not recognized as an up arrow ({up_arrow}),left_arrow ({left_arrow}), up_left_arrow ({up_left_arrow}), or a stop ({stop}).")
        # print(traceback_array,-row,-col,traceback_array[-row,-col])
        print(aligned_seq1)
        print(alignment_indicator)
        print(aligned_seq2)

    return aligned_seq1, aligned_seq2


traceback_alignment(traceback_array, seq1, seq2)
