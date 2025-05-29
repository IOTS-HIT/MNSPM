import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class MNSPM_P {

    static int K = 80000; // The max sequence number of sequence database
    static int maxlen = 10; // length constraint
    static int mingap, maxgap; // gap constraint
    static int minsup; // minimum support
    static int Maxseqlen = 1; // the maximum length of sequence in the database

    static  class Sequence {
        int id;
        int[] S;
        int length;

        Sequence(int id, int[] S, int length) {
            this.id = id;
            this.S = S;
            this.length = length;
        }
    }

    static class Pattern {
        int[] pattern;
        int length;

        Pattern(int[] pattern, int length) {
            this.pattern = pattern;
            this.length = length;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Pattern that = (Pattern) o;
            return Arrays.equals(pattern, that.pattern);
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(pattern);
        }

        @Override
        public String toString() {
            return Arrays.toString(pattern);
        }
    }

    static class Candidate_table {
        HashMap<Integer,List<Integer>> candidate = new HashMap<>(); // store the candidate table

        Candidate_table() {
            this.candidate = new HashMap<>();
        }
    }

    static class HBTable {
        TreeMap<Integer, List<List<Integer>>> HBT;

        HBTable() {
            this.HBT = new TreeMap<>();
        }
    }

    static class HBTableSet {
        HBTable[] HBTS;

        HBTableSet(int length) {
            this.HBTS = new HBTable[length];
            for (int i = 0; i < length; i++) {
                this.HBTS[i] = new HBTable();
            }
        }
    }

    static int SDB_num = 0; // real number of sequences in the database
    static Sequence[] SDB = new Sequence[K];  // Sequence database
    static int AlphabetSize = 1; // the size of alphabet
    static int Cut3_num = 0; // the number of 3-pattern that are cutted
    static int Cut4_num = 0; // the number of 4-pattern that are cutted
    static Map<Integer, Integer> Alphabet = new LinkedHashMap<>(); // store the alphabet of the sequence database
    static List<Pattern> candidate = new ArrayList<>(); // store candidate patterns
    static Map<Integer, List<List<Integer>>> temp_HBT = new HashMap<>(); // store the HBTableSet of each pattern in each sequence
    static Map<Pattern, HBTableSet> All_HBTS = new HashMap<>(); // store the HBTableSet of all patterns in each sequence
    static List<Integer> position = new ArrayList<>(); // store the position of candidate patterns
    static List<Integer>[] record = new ArrayList[K]; // record the position used already in each position of pattern

    static List<Pattern> ReadFile(Sequence[] SDB, double reminsup, String filename) {
        // Read the sequence database from the file
        List<Pattern> frequent_item = new ArrayList<>();
        try (Scanner fileScanner = new Scanner(new File(filename), "UTF-8")) {
//            System.out.println("File opened successfully.");
            int i = 0, length = 0;
            while (fileScanner.hasNextLine() && i < K) {
                String line = fileScanner.nextLine();
                String[] S = line.split("\\s+");
                int[] intArray = java.util.Arrays.stream(S)
                        .filter(s -> !"-1".equals(s) && !"-2".equals(s) && !s.isEmpty())
                        .mapToInt(Integer::parseInt)
                        .toArray();
                SDB[i] = new Sequence(i, intArray, intArray.length);
                length += intArray.length;
                if(SDB[i].length > Maxseqlen)
                    Maxseqlen = SDB[i].length;
                for (int item : intArray)
                    Alphabet.put(item, Alphabet.getOrDefault(item, 0) + 1);
                i++;
            }
            SDB_num = i; // real number of sequences in the database
            minsup = (int) (reminsup * SDB_num);
            AlphabetSize = Alphabet.size();
//            System.out.println("Average length of sequence database: " + (float)length/SDB_num);
            for (Map.Entry<Integer, Integer> entry : Alphabet.entrySet()) {
                int item = entry.getKey();
                int count = entry.getValue();
                if (count >= minsup) {
                    frequent_item.add(new Pattern(new int[]{item}, 1));
                }
            }
//            System.out.println("Frequent item set: " + frequent_item);
//            System.out.println("Frequent item size: " + frequent_item.size());
        } catch (FileNotFoundException e) {
            System.out.println("Failed to open the file.");
            System.exit(0);
        }
        return frequent_item;
    }

    static void Build_PositionMap(Sequence[] SDB, List<Pattern> frequent_item, List<Integer>[] PositionList, Map<Integer, Integer>[] IndexMap) {

        for (int i = 0; i < SDB_num; i++) {
            Map<Integer, List<Integer>> PositionArray = new HashMap<>();
            for (int j = 0; j < SDB[i].length; j++) {
                int item = SDB[i].S[j];
                PositionArray.computeIfAbsent(item, k -> new ArrayList<>()).add(j);
            }
            int len = 0;
            for (int item : PositionArray.keySet()) {
                int flag = 0;
                for (Pattern p : frequent_item)
                    if (p.pattern[0] == item){
                        flag = 1;
                        break;
                    }
                if(flag == 1 && !PositionArray.get(item).isEmpty()){
                    PositionList[i].addAll(PositionArray.get(item));
                    IndexMap[i].put(item, len);
                    len += PositionArray.get(item).size();
                }
            }
        }
    }

    static void Generate_candidate(List<Pattern> candidate, int length, List<Pattern>[] frequent) {
        // Generate (length+1)—candidate patterns based on length-frequent patterns
        if (length == 1) {
            for (Pattern p : frequent[0]) {
                for (Pattern q : frequent[0]) {
                    Pattern super_pattern = new Pattern(new int[]{p.pattern[0], q.pattern[0]}, 2);
                    candidate.add(super_pattern);
                }
            }
        }
        else {
            for (Pattern p : frequent[length - 1]) {
                for (Pattern q : frequent[length - 1]) {
                    int flag = 0;
                    for (int i = 0; i < length - 1; i++){
                        if (p.pattern[i] != q.pattern[i+1]){
                            flag = 1;
                            break;
                        }
                    }
                    if (flag == 0){
                        int[] pattern = Arrays.copyOf(q.pattern, length+1);
                        pattern[length] = p.pattern[length-1];
                        Pattern super_pattern = new Pattern(pattern, length+1);
                        candidate.add(super_pattern);
                    }
                }
            }
        }
    }

    static void GapMatch(List<Integer>[] PositionList, Map<Integer, Integer>[] IndexMap, Pattern candidate, List<Pattern>[] frequent) {
        // Mining 2-frequent patterns based on 2-candidate patterns and store all HBTables
        int first_item = candidate.pattern[0];
        int second_item = candidate.pattern[1];
        int support = 0;
        HBTableSet HBTS = new HBTableSet(SDB_num);

        for (int i = 0; i < SDB_num; i++) {
//            System.out.println("Processing sequence " + i);
            if(SDB[i].length == 0)  continue;
            List<Integer> first_pos = new ArrayList<>();
            List<Integer> second_pos = new ArrayList<>();
            if (!IndexMap[i].containsKey(first_item) || !IndexMap[i].containsKey(second_item))
                continue;
            // get the positions of first and second item in the sequence
            int found = 0;
            for (Integer key : IndexMap[i].keySet()) {
                if (found == 1) {
                    first_pos.addAll(PositionList[i].subList(IndexMap[i].get(first_item), IndexMap[i].get(key)));
                    found++;
                    break;
                }
                if (key == first_item)
                    found = 1;
            }
            if(found == 1)
                first_pos.addAll(PositionList[i].subList(IndexMap[i].get(first_item), PositionList[i].size()));

            found = 0;
            for (Integer key : IndexMap[i].keySet()) {
                if (found == 1) {
                    second_pos.addAll(PositionList[i].subList(IndexMap[i].get(second_item), IndexMap[i].get(key)));
                    found++;
                    break;
                }
                if (key == second_item)
                    found = 1;
            }
            if(found == 1)
                second_pos.addAll(PositionList[i].subList(IndexMap[i].get(second_item), PositionList[i].size()));

            // build the candidate table
            Candidate_table CT = new Candidate_table();
            for (Integer pos : first_pos) {
                for (int j = mingap; j <= maxgap; j++) {
                    int cand = pos + j + 1;
                    if (CT.candidate.containsKey(cand))
                        CT.candidate.get(cand).add(pos);
                    else
                        CT.candidate.put(cand, new ArrayList<>(Arrays.asList(pos)));
                }
            }

            // stroe the matched positions into HBT
            for (int cand : second_pos) {
                if (!CT.candidate.containsKey(cand))
                    continue;
                for (int pos : CT.candidate.get(cand)) {
                    if (!HBTS.HBTS[i].HBT.containsKey(pos)) {
                        List<Integer> list = new ArrayList<>();
                        list.add(cand);
                        List<List<Integer>> list_list = new ArrayList<>();
                        list_list.add(list);
                        HBTS.HBTS[i].HBT.put(pos, list_list);
                    } else {
                        List<Integer> list = new ArrayList<>();
                        list.add(cand);
                        HBTS.HBTS[i].HBT.get(pos).add(list);
                    }
                }
            }

            // get all combinations of positions and filter out the non-overlap combinations
            for (Integer pos : HBTS.HBTS[i].HBT.keySet()) {
                if (HBTS.HBTS[i].HBT.get(pos).isEmpty()) {
                    HBTS.HBTS[i].HBT.remove(pos);
                    continue;
                }
                for (List<Integer> occ : HBTS.HBTS[i].HBT.get(pos)) {
                    if (!(record[0].contains(pos) || record[1].contains(occ.get(0)))) {
                        support++;
                        record[0].add(pos);
                        record[1].add(occ.get(0));
                        break;
                    }
                }
            }
            position.clear();
            record[0].clear();
            record[1].clear();
        }
        if (support >= minsup) {
            frequent[1].add(candidate);
            All_HBTS.put(candidate, HBTS);
        }
    }

    static void Merge_HBTable(HBTable HBT1, HBTable HBT2, Map<Integer, List<List<Integer>>> temp_HBT) {
        // Merge two HBTables and store the result in temp_HBT
        for (int pos : HBT1.HBT.keySet()) {
            for (List<Integer> left_occ : HBT1.HBT.get(pos)) {
                if (!HBT2.HBT.containsKey(left_occ.get(left_occ.size() - 1)))
                    continue;
                for (List<Integer> right_occ : HBT2.HBT.get(left_occ.get(left_occ.size() - 1))) {
                    List<Integer> occ = new ArrayList<>();
                    occ.add(left_occ.get(0));
                    occ.add(right_occ.get(0));
                    temp_HBT.get(pos).add(occ);
                }
            }
        }
    }

    static int Calculate_support(HBTableSet HBTS1, HBTableSet HBTS2, int length) {
        // Calculate the support of a HBTableSet
        int support = 0;
        for (int i = 0; i < SDB_num; i++) {
            outerLoop:
            for (int pos : HBTS1.HBTS[i].HBT.keySet()) {
                for (List<Integer> left_occ : HBTS1.HBTS[i].HBT.get(pos)) {
                    int common_pos = left_occ.get(left_occ.size() - 1);
                    if (!HBTS2.HBTS[i].HBT.containsKey(common_pos))
                        continue;
                    for (List<Integer> right_occ : HBTS2.HBTS[i].HBT.get(common_pos)) {
                        int flag = 0;
                        if (record[0].contains(pos) || record[length-1].contains(right_occ.get(0)))
                            flag = 1;
                        else
                            for (int l = 1; l < length - 1; l++)
                                if (record[l].contains(left_occ.get(l - 1))) {
                                    flag = 1;
                                    break;
                                }
                        if (flag == 0) {
                            support++;
                            record[0].add(pos);
                            for (int l = 1; l < length - 1; l++)
                                record[l].add(left_occ.get(l - 1));
                            record[length-1].add(right_occ.get(0));
                            continue outerLoop;
                        }
                    }
                }
            }
            for (int k = 0; k < length; k++)
                record[k].clear();
            if (support >= minsup)
                break;
        }
        return support;
    }

    static int Calculate_cutnum(Pattern candidate, int length) {
        int cut_num = 0;
        Pattern[] substr = new Pattern[length - 1];
        for (int i = 0; i < length - 1; i++) {
            substr[i] = new Pattern(new int[]{candidate.pattern[i], candidate.pattern[i + 1]}, 2);
        }

        if (length == 3) {
            for (int i = 0; i < SDB_num; i++) {
                List<Integer> pos_list = new ArrayList<>();
                Pattern substr0 = new Pattern(substr[0].pattern, substr[0].length);
                Pattern substr1 = new Pattern(substr[1].pattern, substr[1].length);

                if (All_HBTS.containsKey(substr0) && All_HBTS.containsKey(substr1)) {
                    for (int pos : All_HBTS.get(substr0).HBTS[i].HBT.keySet()) {
                        for (List<Integer> left_occ : All_HBTS.get(substr0).HBTS[i].HBT.get(pos)) {
                            if (!pos_list.contains(left_occ.get(left_occ.size() - 1))) {
                                pos_list.add(left_occ.get(left_occ.size() - 1));
                            }
                        }
                    }
                    List<Integer> pos_list2 = new ArrayList<>(All_HBTS.get(substr1).HBTS[i].HBT.keySet());
                    pos_list.retainAll(pos_list2);
                    cut_num += pos_list.size();
                    if (cut_num >= minsup) {
                        break;
                    }
                }
            }
        } else if (length == 4) {
            for (int i = 0; i < SDB_num; i++) {
                List<Integer> pos_list = new ArrayList<>();
                Pattern substr0 = new Pattern(substr[0].pattern, substr[0].length);
                Pattern substr1 = new Pattern(substr[1].pattern, substr[1].length);
                Pattern substr2 = new Pattern(substr[2].pattern, substr[2].length);

                if (All_HBTS.containsKey(substr0) && All_HBTS.containsKey(substr1) && All_HBTS.containsKey(substr2)) {
                    for (int pos : All_HBTS.get(substr0).HBTS[i].HBT.keySet()) {
                        for (List<Integer> left_occ : All_HBTS.get(substr0).HBTS[i].HBT.get(pos)) {
                            if (!All_HBTS.get(substr1).HBTS[i].HBT.containsKey(left_occ.get(left_occ.size() - 1))) {
                                continue;
                            }
                            for (List<Integer> right_occ : All_HBTS.get(substr1).HBTS[i].HBT.get(left_occ.get(left_occ.size() - 1))) {
                                if (!pos_list.contains(right_occ.get(right_occ.size() - 1))) {
                                    pos_list.add(right_occ.get(right_occ.size() - 1));
                                }
                            }
                        }
                    }
                    List<Integer> pos_list2 = new ArrayList<>(All_HBTS.get(substr2).HBTS[i].HBT.keySet());
                    pos_list.retainAll(pos_list2);
                    cut_num += pos_list.size();
                    if (cut_num >= minsup) {
                        break;
                    }
                }
            }
        } else {
            Map<Integer, List<List<Integer>>> HBT = new HashMap<>();
            for (int i = 0; i < Maxseqlen; i++) {
                HBT.put(i, new ArrayList<>());
            }
            for (int i = 0; i < SDB_num; i++) {
                Pattern substr0 = new Pattern(substr[0].pattern, substr[0].length);
                Pattern substr1 = new Pattern(substr[1].pattern, substr[1].length);

                if (All_HBTS.containsKey(substr0) && All_HBTS.containsKey(substr1)) {
                    Merge_HBTable(All_HBTS.get(substr0).HBTS[i], All_HBTS.get(substr1).HBTS[i], HBT);
                    for (int l = 2; l < length - 2; l++) {
                        Pattern substr_l = new Pattern(substr[l].pattern, substr[l].length);
                        for (int pos : All_HBTS.get(substr0).HBTS[i].HBT.keySet()) {
                            List<List<Integer>> temp_list = new ArrayList<>();
                            for (List<Integer> occ : HBT.get(pos)) {
                                List<Integer> temp = new ArrayList<>(occ);
                                temp_list.add(temp);
                            }
                            HBT.get(pos).clear();
                            for (List<Integer> left_occ : temp_list) {
                                if (!All_HBTS.get(substr_l).HBTS[i].HBT.containsKey(left_occ.get(left_occ.size() - 1))) {
                                    continue;
                                }
                                for (List<Integer> right_occ : All_HBTS.get(substr_l).HBTS[i].HBT.get(left_occ.get(left_occ.size() - 1))) {
                                    List<Integer> occ = new ArrayList<>(left_occ);
                                    occ.add(right_occ.get(0));
                                    HBT.get(pos).add(occ);
                                }
                            }
                        }
                    }
                    List<Integer> pos_list = new ArrayList<>();
                    for (int pos : All_HBTS.get(substr0).HBTS[i].HBT.keySet()) {
                        for (List<Integer> left_occ : HBT.get(pos)) {
                            if (!pos_list.contains(left_occ.get(left_occ.size() - 1))) {
                                pos_list.add(left_occ.get(left_occ.size() - 1));
                            }
                        }
                    }
                    Pattern substr_last = new Pattern(substr[length - 2].pattern, substr[length - 2].length);
                    List<Integer> pos_list2 = new ArrayList<>(All_HBTS.get(substr_last).HBTS[i].HBT.keySet());
                    pos_list.retainAll(pos_list2);
                    cut_num += pos_list.size();
                    for (int j = 0; j < length; j++) {
                        HBT.get(j).clear();
                    }
                    if (cut_num >= minsup) {
                        break;
                    }
                }
            }
        }
        return cut_num;
    }

    static void Mining_patterns_withcut(Sequence[] SDB, Pattern candidate, List<Pattern>[] frequent) {
        int length = candidate.length;
        int support = 0, Length;
        Pattern[] substr = new Pattern[length - 1];
        for (int i = 0; i < length - 1; i++) {
            substr[i] = new Pattern(new int[]{candidate.pattern[i], candidate.pattern[i + 1]}, 2);
        }

        // 挖掘3-频繁模式
        if (length == 3) {
            int cut = Calculate_cutnum(candidate, 3);
            if (cut >= minsup) {
                Pattern substr0 = new Pattern(substr[0].pattern, substr[0].length);
                Pattern substr1 = new Pattern(substr[1].pattern, substr[1].length);
                support = Calculate_support(All_HBTS.get(substr0), All_HBTS.get(substr1), 3);
                if (support >= minsup) {
                    frequent[length - 1].add(candidate);
                }
            } else {
                Cut3_num++;
            }
        }

        // 挖掘4-频繁模式
        if (length == 4) {
            int cut = Calculate_cutnum(candidate, 4);
            if (cut >= minsup) {
                for (int i = 0; i < SDB_num; i++) {
                    Length = SDB[i].length;
                    Pattern substr0 = new Pattern(substr[0].pattern, substr[0].length);
                    Pattern substr1 = new Pattern(substr[1].pattern, substr[1].length);
                    Pattern substr2 = new Pattern(substr[2].pattern, substr[2].length);

                    if (All_HBTS.containsKey(substr0) && All_HBTS.containsKey(substr1) && All_HBTS.containsKey(substr2)) {
                        Merge_HBTable(All_HBTS.get(substr0).HBTS[i], All_HBTS.get(substr1).HBTS[i], temp_HBT);
                        outerLoop:
                        for (int pos : All_HBTS.get(substr0).HBTS[i].HBT.keySet()) {
                            for (List<Integer> left_occ : temp_HBT.get(pos)) {
                                if (!All_HBTS.get(substr2).HBTS[i].HBT.containsKey(left_occ.get(left_occ.size() - 1)))
                                    continue;
                                for (List<Integer> right_occ : All_HBTS.get(substr2).HBTS[i].HBT.get(left_occ.get(left_occ.size() - 1))) {
                                    int flag = 0;
                                    if (record[0].contains(pos) || record[length - 1].contains(right_occ.get(0)))
                                        flag = 1;
                                    else {
                                        for (int l = 1; l < length - 1; l++) {
                                            if (record[l].contains(left_occ.get(l - 1))) {
                                                flag = 1;
                                                break;
                                            }
                                        }
                                    }
                                    if (flag == 0) {
                                        support++;
                                        record[0].add(pos);
                                        for (int l = 1; l < length - 1; l++)
                                            record[l].add(left_occ.get(l - 1));
                                        record[length - 1].add(right_occ.get(0));
                                        continue outerLoop;
                                    }
                                }
                            }
                        }

                        for (int j = 0; j < Length; j++)
                            temp_HBT.get(j).clear();
                        for (int j = 0; j < length; j++)
                            record[j].clear();
                        if (support >= minsup) {
                            frequent[length - 1].add(candidate);
                            break;
                        }
                    }
                }
            } else {
                Cut4_num++;
            }
        }

        // 针对5-频繁模式或更长的模式，进行更复杂的挖掘
        if (length >= 5) {
            int cut = Calculate_cutnum(candidate, length);
            if (cut >= minsup) {
                for (int i = 0; i < SDB_num; i++) {
                    Length = SDB[i].length;
                    Pattern substr0 = new Pattern(substr[0].pattern, substr[0].length);
                    Pattern substr1 = new Pattern(substr[1].pattern, substr[1].length);

                    if (All_HBTS.containsKey(substr0) && All_HBTS.containsKey(substr1)) {
                        Merge_HBTable(All_HBTS.get(substr0).HBTS[i], All_HBTS.get(substr1).HBTS[i], temp_HBT);
                        for (int l = 2; l < length - 2; l++) {
                            Pattern substr_l = new Pattern(substr[l].pattern, substr[l].length);
                            for (int pos : All_HBTS.get(substr0).HBTS[i].HBT.keySet()) {
                                List<List<Integer>> temp_list = new ArrayList<>();
                                for (List<Integer> occ : temp_HBT.get(pos)) {
                                    List<Integer> temp = new ArrayList<>(occ);
                                    temp_list.add(temp);
                                }
                                temp_HBT.get(pos).clear();
                                for (List<Integer> left_occ : temp_list) {
                                    if (!All_HBTS.get(substr_l).HBTS[i].HBT.containsKey(left_occ.get(left_occ.size() - 1)))
                                        continue;
                                    for (List<Integer> right_occ : All_HBTS.get(substr_l).HBTS[i].HBT.get(left_occ.get(left_occ.size() - 1))) {
                                        List<Integer> occ = new ArrayList<>(left_occ);
                                        occ.add(right_occ.get(0));
                                        temp_HBT.get(pos).add(occ);
                                    }
                                }
                            }
                        }
                        outerLoop:
                        for (int pos : All_HBTS.get(substr0).HBTS[i].HBT.keySet()) {
                            for (List<Integer> left_occ : temp_HBT.get(pos)) {
                                if (!All_HBTS.get(substr[length - 2]).HBTS[i].HBT.containsKey(left_occ.get(left_occ.size() - 1)))
                                    continue;
                                for (List<Integer> right_occ : All_HBTS.get(substr[length - 2]).HBTS[i].HBT.get(left_occ.get(left_occ.size() - 1))) {
                                    int flag = 0;
                                    if (record[0].contains(pos) || record[length - 1].contains(right_occ.get(0)))
                                        flag = 1;
                                    else {
                                        for (int l = 1; l < length - 1; l++) {
                                            if (record[l].contains(left_occ.get(l - 1))) {
                                                flag = 1;
                                                break;
                                            }
                                        }
                                    }
                                    if (flag == 0) {
                                        support++;
                                        record[0].add(pos);
                                        for (int l = 1; l < length - 1; l++)
                                            record[l].add(left_occ.get(l - 1));
                                        record[length - 1].add(right_occ.get(0));
                                        continue outerLoop;
                                    }
                                }
                            }
                        }

                        for (int j = 0; j < Length; j++)
                            temp_HBT.get(j).clear();
                        for (int j = 0; j < length; j++)
                            record[j].clear();
                        if (support >= minsup) {
                            frequent[length - 1].add(candidate);
                            break;
                        }
                    }
                }
            }
        }
    }


    public static void main(String[] args) {

        // 为方便测试，输入参数组合整理如下，按索引选择即可
        String[] FileDir = new String[6];

        FileDir[0] = "DataSet/Book1.txt";
        FileDir[1] = "DataSet/Book2.txt";
        FileDir[2] = "DataSet/BMS1.txt";
        FileDir[3] = "DataSet/BMS2.txt";
        FileDir[4] = "DataSet/Sign.txt";
        FileDir[5] = "DataSet/Bike.txt";

        System.out.println("Welcome to the MNSPM-P");
//        Scanner scanner = new Scanner(System.in);
//        int number = scanner.nextInt();
//        mingap = scanner.nextInt();
//        maxgap = scanner.nextInt();
////        minsup = scanner.nextInt();
//        double reminsup = scanner.nextDouble();
//        scanner.close();
        int number = Integer.parseInt(args[0]);
        mingap = Integer.parseInt(args[1]);
        maxgap = Integer.parseInt(args[2]);
        double reminsup = Double.parseDouble(args[3]);

        // read the data set and set the parameters
        List<Pattern> frequent_item = ReadFile(SDB, reminsup, FileDir[number]);
//        System.out.println("The minimum support is: " + minsup);

        // initialize the data structures
        List<Integer>[] PositionList = new ArrayList[SDB_num];
        for (int i = 0; i < SDB_num; i++)
            PositionList[i] = new ArrayList<>();
        Map<Integer,Integer>[] IndexList = new HashMap[SDB_num];
        for (int i = 0; i < SDB_num; i++)
            IndexList[i] = new LinkedHashMap<>();
        // store frequent patterns
        List<Pattern>[] frequent = new ArrayList[maxlen];
        frequent[0] = frequent_item;
        for (int i = 1; i < maxlen; i++)
            frequent[i] = new ArrayList<>();
        for (int i = 0; i < maxlen; i++)
            record[i] = new ArrayList<>();
        for (int i = 0; i < Maxseqlen; i++)
            temp_HBT.put(i, new ArrayList<>());

        // monitor the memory usage
        MemoryUsageMonitor monitor = new MemoryUsageMonitor();
        monitor.start();
        // record the time usage
        int begintime = (int) System.currentTimeMillis();

        // build the PoistionMap
        Build_PositionMap(SDB, frequent_item, PositionList, IndexList);
        // mining frequent 2-patterns
        int candidate_num = 0;
        Generate_candidate(candidate, 1, frequent);
//        System.out.println(2 + "-Candidate Patterns: " + candidate.size());
        candidate_num += candidate.size();
        for (Pattern cand : candidate) {
            GapMatch(PositionList, IndexList, cand, frequent);
        }
        candidate.clear();

        // mining frequent 3-patterns and longer patterns
        int length = 2;
        int cand3 = 0;
        int cand4 = 0;
        while (length < maxlen && !frequent[length - 1].isEmpty()) {
            Generate_candidate(candidate, length, frequent);
            candidate_num += candidate.size();
            // System.out.println(length+1 + "-Candidate Patterns: " + candidate.size());
            length++;
            if (length == 3)
                cand3 = candidate.size();
            if (length == 4)
                cand4 = candidate.size();

            for (Pattern cand : candidate) {
                Mining_patterns_withcut(SDB, cand, frequent);
            }
            candidate.clear();
        }

        // get the time usage
        int endtime = (int) System.currentTimeMillis();
        // get the memory usage
        double maxMemoryUsage = monitor.getMaxMemoryUsage();
        monitor.stopMonitoring();

        // output the frequent patterns
        System.out.println(FileDir[number]);
        int pattern_num = 0;
        for (int l = 0; l < maxlen; l++) {
            if (frequent[l].isEmpty())
                break;
//            System.out.println(l + "-frequent patterns: " + frequent[l].size());
            pattern_num += frequent[l].size();
        }

//        System.out.println("3-candidate pattern num: " + cand3);
//        System.out.println("3-frequent pattern num: " + frequent[2].size());
//        System.out.println("3-infrequent pattern num: " + (cand3 - frequent[2].size()));
//        System.out.println("3-candidiate-cutted number: " + Cut3_num);
//        System.out.println("3-模式剪枝率: " + (double)Cut3_num/(cand3 - frequent[2].size()));
//
//        System.out.println("4-candidate pattern num: " + cand4);
//        System.out.println("4-frequent pattern num: " + frequent[3].size());
//        System.out.println("4-infrequent pattern num: " + (cand4 - frequent[3].size()));
//        System.out.println("4-candidiate-cutted number: " + Cut4_num);
//        System.out.println("4-模式剪枝率: " + (double)Cut4_num/(cand4 - frequent[3].size()));
//
        System.out.println("The number of candidate patterns: " + candidate_num);
        System.out.println("The number of frequent patterns: " + pattern_num);
        System.out.println("The number of infrequent patterns:" + (candidate_num - pattern_num));

        // output the result
        DecimalFormat df = new DecimalFormat("#.##");
        System.out.println("Time used: " + df.format((endtime - begintime) / 1000.0) + "s");
        System.out.println("Memory used: " + df.format(maxMemoryUsage) + "MB");

    }
}
