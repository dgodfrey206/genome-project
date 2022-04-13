#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

class Blast_DB {
 public:
  Blast_DB(std::string genome)
      : genome_(std::move(genome)), seed_pos(10'000'000) {
    
  }

  ~Blast_DB() = default;

  auto& table() { return seed_pos; }

  struct data {
    std::string polymer;
    std::size_t query_index;
    std::size_t genome_index;
  };
  static const int WORD_SIZE = 11;

 private:
  UnorderedMapPool seed_pos;
  std::vector<data> stk;
  std::string genome_;

  static const int ORIGINAL_SIZE = 50;

 public:
  static auto traceback_alignment(
      std::vector<std::vector<std::string>> traceback_array,
      std::string seq1,
      std::string seq2,
      std::string up_arrow = "↑",
      std::string left_arrow = "←",
      std::string up_left_arrow = "🡔",
      std::string stop = "-") {
    /*Align seq1 and seq2 using the traceback matrix and return as two strings

    traceback_array -- a numpy array with arrow characters indicating the
    direction from
    which the best path to a given alignment position originated

    seq1 - a sequence represented as a string
    seq2 - a sequence represented as a string
    up_arrow - the unicode used for the up arrows (there are several arrow
    symbols in Unicode)
    left_arrow - the unicode used for the left arrows
    up_left_arrow - the unicode used for the diagonal arrows
    stop - the symbol used in the upper left to indicate the end of the
    alignment
    */

    auto n_rows = seq1.size() + 1;     // need an extra row up top
    auto n_columns = seq2.size() + 1;  // need an extra row up top

    auto row = seq1.size();
    auto col = seq2.size();
    auto arrow = traceback_array[row][col];
    std::string aligned_seq1 = "";
    std::string aligned_seq2 = "";
    std::string alignment_indicator = "";
    while (arrow != "-") {
      arrow = traceback_array[row][col];

      if (arrow == up_arrow) {
        // We want to add the new indel onto the left
        // side of the growing aligned sequence
        aligned_seq2 = "-" + aligned_seq2;
        aligned_seq1 = seq1[row - 1] + aligned_seq1;
        alignment_indicator = " " + alignment_indicator;
        row -= 1;
      } else if (arrow == up_left_arrow) {
        // Note that we look up the row-1 and col-1 indexes
        // because there is an extra "-" character at the
        // start of each sequence
        auto seq1_character = seq1[row - 1];
        auto seq2_character = seq2[col - 1];
        aligned_seq1 = seq1[row - 1] + aligned_seq1;
        aligned_seq2 = seq2[col - 1] + aligned_seq2;
        if (seq1_character == seq2_character)
          alignment_indicator = "|" + alignment_indicator;
        else
          alignment_indicator = " " + alignment_indicator;
        row -= 1;
        col -= 1;
      } else if (arrow == left_arrow) {
        aligned_seq1 = "-" + aligned_seq1;
        aligned_seq2 = seq2[col - 1] + aligned_seq2;
        alignment_indicator = " " + alignment_indicator;
        col -= 1;
      } else if (arrow == stop)
        break;
      else {
        auto s="Traceback array entry at " + std::to_string(row) + "," + std::to_string(col) + ": " + arrow + " is not recognized "
            "as an up arrow " + (up_arrow) + ",left_arrow " + (left_arrow) + ", "
            "up_left_arrow " + (up_left_arrow) + ", or a stop " + (stop) + ".";
        std::cerr << s << std::endl;
        throw std::invalid_argument("Error");
      }
    }
    
    return std::make_pair(aligned_seq1, aligned_seq2);
  }

  static auto query(std::string const& seq1, std::string const& seq2, vector<vector<int>>* s = 0, vector<vector<std::string>>* t = 0) {
    // build an array of zeroes
    // seq1 = "GCTGATTC"
    // seq2 = "GATCTGATTA"
    auto n_rows = 1 + seq1.size();
    auto n_columns = 1 + seq2.size();

    std::vector<std::vector<int>> scoring_array(n_rows, std::vector<int>(n_columns, 0));
    std::vector<std::vector<std::string>> traceback_array(n_rows, std::vector<std::string>(n_columns, "-"));

    int count = 0;
    for (int i = 0; i < n_rows; i++) {
      for (int j = 0; j < n_columns; j++) {
        scoring_array[i][j] = count;
        count += 1;
      }
    }

    // pretty_table_from_array(scoring_array, row_labels, column_labels)

    std::string up_arrow = "↑";
    std::string right_arrow = "→";
    std::string down_arrow = "↓";
    std::string left_arrow = "←";
    std::string down_right_arrow = "🡖";
    std::string up_left_arrow = "🡔";
 
    // build an array of zeroes
    n_rows = seq1.size() + 1;     // need an extra row up top
    n_columns = seq2.size() + 1;  // need an extra column on the left

    scoring_array.clear();
    scoring_array.resize(n_rows, std::vector<int>(n_columns, 0));
    traceback_array.clear();
    traceback_array.resize(n_rows, std::vector<std::string>(n_columns, "-"));

    std::string arrow = "-";
    int gap_penalty = -1;
    int match_bonus = 2;
    int mismatch_penalty = -1;

    int score = 0;
    int previous_score = 0;

    int cell_to_the_left = 0;
    int from_left_score = 0;
    int above_cell = 0;
    int from_above_score = 0;
    int diagonal_left_cell = 0;
    int diagonal_left_cell_score = 0;
    // iterate over columns first because we want to do
    // all the columns for row 1 before row 2
    int row = 0,col = 0;
    for (row = 0; row < n_rows; row++) {
      for (col = 0; col < n_columns; col++) {
        if (row == 0 && col == 0) {
          score = 0;
          arrow = "-";
        } else if (row == 0) {
          previous_score = scoring_array.at(row).at(col - 1);
          
          score = previous_score + gap_penalty;
          arrow = left_arrow;
        } else if (col == 0) {
          // We're on the first column but not in the first row
          previous_score = scoring_array.at(row - 1).at(col);
          score = previous_score + gap_penalty;
          arrow = up_arrow;
          
        } else {
          cell_to_the_left = scoring_array.at(row).at(col - 1);
          from_left_score = cell_to_the_left + gap_penalty;

          above_cell = scoring_array.at(row - 1).at(col);
          from_above_score = above_cell + gap_penalty;

          diagonal_left_cell = scoring_array.at(row - 1).at(col - 1);

          if (seq1.at(row - 1) == seq2.at(col - 1))
            diagonal_left_cell_score = diagonal_left_cell + match_bonus;
          else
            diagonal_left_cell_score = diagonal_left_cell + mismatch_penalty;

          score = std::max(
              {from_left_score, from_above_score, diagonal_left_cell_score});
              
          if (score == from_left_score)
            arrow = left_arrow;
          else if (score == from_above_score)
            arrow = up_arrow;
          else if (score == diagonal_left_cell_score)
            arrow = up_left_arrow;
        }
        traceback_array.at(row).at(col) = arrow;
        scoring_array.at(row).at(col) = score;
      }
    }

    if (s) *s = scoring_array;
    if (t) *t = traceback_array;

    return std::make_pair(scoring_array[row-1][col-1], traceback_alignment(traceback_array, seq1, seq2));
    //return traceback_alignment(traceback_array, scoring_array, seq1, seq2);
  }

  void store_polymers() {
    UnorderedMapPool found;
    for (int i = 0; i < genome_.size(); i++) {
      std::string word = genome_.substr(i, WORD_SIZE);
      if (word.size() == WORD_SIZE && found[word] == 0) {
        found[word] = 1;
        seed_pos[word] = i - WORD_SIZE;
      }
    }
  }
};