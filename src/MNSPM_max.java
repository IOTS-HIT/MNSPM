// This is the main program of the VLG-PM algorithm made by MXR
import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class MNSPM_max {

    static int K = 1200; // The max sequence number of sequence database
    static int maxlen = 10; // length constraint
    static int mingap, maxgap; // gap constraint
    static int minsup; // minimum support
    static int Maxseqlen = 1; // the maximum length of sequence in the database

    static class Sequence {
        int id; // sequence id
        char[] S; // sequence
        int length; // length of sequence

        Sequence(int id, char[] S, int length) {
            this.id = id;
            this.S = S;
            this.length = length;
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
    static int AlphabetSize = 26; // the size of alphabet
    static List<String> candidate = new ArrayList<>(); // store candidate patterns
    static List<Integer> position = new ArrayList<>(); // store the position of candidate patterns
    static List<Integer>[] record = new ArrayList[K]; // record the position used already in each position of pattern
    static Map<Integer, List<List<Integer>>> temp_HBT = new HashMap<>(); // store the HBTableSet of each pattern in each sequence
    static Map<String, HBTableSet> All_HBTS = new HashMap<>(); // store the HBTableSet of all patterns in each sequence

    static void ReadFile(Sequence[] SDB, String filename) {
        // Read the sequence database from the file
        try (Scanner fileScanner = new Scanner(new File(filename), "UTF-8")) {
//            System.out.println("File opened successfully.");
            int i = 0;
            int length = 0;
            while (fileScanner.hasNextLine() && i < K) {
                // read every single line of the file and store it in the sequence database
                String line = fileScanner.nextLine();
                SDB[i] = new Sequence(i, line.toCharArray(), line.length());
                length += line.length();
                if(SDB[i].length > Maxseqlen)
                    Maxseqlen = SDB[i].length;
                i++;
            }
            SDB_num = i; // real number of sequences in the database
//            System.out.println("Average length of sequence database: " + (float)length/SDB_num);
        } catch (FileNotFoundException e) {
            System.out.println("Failed to open the file.");
            System.exit(0);
        }
    }

    static void Build_PositionMap(Sequence[] SDB, List<Integer>[] PositionList, Map<Integer, Integer>[] IndexMap) {

        List<Integer>[] Temp_PositionArray = new ArrayList[AlphabetSize];
        for (int i = 0; i < AlphabetSize; i++)
            Temp_PositionArray[i] = new ArrayList<>();

        for (int i = 0; i < SDB_num; i++) {
            for (int j = 0; j < SDB[i].length; j++) {
                if (SDB[i].S[j] >= 'a' && SDB[i].S[j] <= 'z')
                    Temp_PositionArray[SDB[i].S[j] - 'a'].add(j);
//                if (SDB[i].S[j] >= 'A' && SDB[i].S[j] <= 'Z')
//                    Temp_PositionArray[SDB[i].S[j] - 'A' + 26].add(j);
            }
            int len = 0;
            for (int j = 0; j < AlphabetSize; j++) {
                if(!Temp_PositionArray[j].isEmpty()) {
                    PositionList[i].addAll(Temp_PositionArray[j]);
                    IndexMap[i].put(j, len);
                    len += Temp_PositionArray[j].size();
                    Temp_PositionArray[j].clear();
                }
                else
                    IndexMap[i].put(j, -1);
            }
        }
    }

    static void Mining_frequent_item(List<Integer>[] PositionList, Map<Integer, Integer>[] IndexMap, List<String>[] frequent) {

        for (int j = 0; j < AlphabetSize; j++) {
            int support = 0;
            for (int i = 0; i < SDB_num; i++) {
                if (IndexMap[i].get(j) != -1) {
                    int temp = j+1;
                    while (temp < AlphabetSize && IndexMap[i].get(temp) == -1)
                        temp++;
                    if (temp < AlphabetSize) {
                        support += IndexMap[i].get(temp) - IndexMap[i].get(j);
                    }
                    else
                        support += PositionList[i].size() - IndexMap[i].get(j);
                    if (support >= minsup) {
                        frequent[0].add(Character.toString((char)(j + 'a')));
                        break;
                    }
                }
            }
        }
    }

    static void Prun_PositionMap(List<Integer>[] PositionList, Map<Integer, Integer>[] IndexMap, List<String>[] frequent) {

        for (int i = 0; i < SDB_num; i++) {
            int len = 0;
            int temp_len = 0;
            for (int j = 0; j < AlphabetSize; j++) {
                if (IndexMap[i].get(j) != -1) {
                    int temp = j+1;
                    while (temp < AlphabetSize && IndexMap[i].get(temp) == -1)
                        temp++;
                    // frequent item -> update index map
                    if (frequent[0].contains(Character.toString((char)(j + 'a')))) {
                        if (temp < AlphabetSize)
                            temp_len = IndexMap[i].get(temp) - IndexMap[i].get(j);
                        IndexMap[i].put(j, len);
                        if (temp < AlphabetSize)
                            len += temp_len;
                    }
                    // infrequent item -> remove from position list and index map
                    else {
                        if (temp < AlphabetSize) {
                            temp_len = IndexMap[i].get(temp) - IndexMap[i].get(j);
                            int end = Math.min(len + temp_len, PositionList[i].size() - 1);
                            PositionList[i].subList(len, end).clear();
                            IndexMap[i].put(j, -1);
                        }
                        else {
                            PositionList[i].subList(len, PositionList[i].size()).clear();
                            IndexMap[i].put(j, -1);
                        }
                    }
                }
            }
        }
    }

    static void Generate_candidate(List<String> candidate,int length, List<String>[] frequent) {
        // Generate (length+1)—candidate patterns based on length-frequent patterns
        if (length == 1) {
            for (int i = 0; i < frequent[0].size(); i++) {
                for (int j = 0; j < frequent[0].size(); j++) {
                    String super_pattern = frequent[0].get(i) + frequent[0].get(j);
                    candidate.add(super_pattern);
                }
            }
        }
        else {
            for (int i = 0; i < frequent[length - 1].size(); i++) {
                for (int j = 0; j < frequent[length - 1].size(); j++) {
                    String suffix = frequent[length - 1].get(i).substring(1);
                    String prefix = frequent[length - 1].get(j).substring(0, length - 1);
                    // pattern connection strategy
                    if (prefix.equals(suffix)) {
                        String super_pattern = frequent[length - 1].get(i).charAt(0) + suffix + frequent[length - 1].get(j).charAt(length - 1);
                        candidate.add(super_pattern);
                    }
                }
            }
        }
    }

    static void GapMatch(List<Integer>[] PositionList, Map<Integer, Integer>[] IndexMap, String candidate, List<String>[] frequent) {
        // Mining 2-frequent patterns based on 2-candidate patterns and store all HBTables
        char first_item = candidate.charAt(0);
        char second_item = candidate.charAt(1);
        int index1 = first_item - 'a';
        int index2 = second_item - 'a';
        int support = 0;
        HBTableSet HBTS = new HBTableSet(SDB_num);

        for (int i = 0; i < SDB_num; i++) {
            // find the first item's position in PositionMap
            List<Integer> first_pos = new ArrayList<>();
            List<Integer> second_pos = new ArrayList<>();
            if (IndexMap[i].get(index1) == -1 || IndexMap[i].get(index2) == -1)
                continue;
            // get the positions of first and second item in the sequence
            int temp1 = index1 + 1;
            while (temp1 < AlphabetSize && IndexMap[i].get(temp1) == -1)
                temp1++;
            if (temp1 < AlphabetSize)
                first_pos.addAll(PositionList[i].subList(IndexMap[i].get(index1), IndexMap[i].get(temp1)));
            else
                first_pos.addAll(PositionList[i].subList(IndexMap[i].get(index1), PositionList[i].size()));
            int temp2 = index2 + 1;
            while (temp2 < AlphabetSize && IndexMap[i].get(temp2) == -1)
                temp2++;
            if (temp2 < AlphabetSize)
                second_pos.addAll(PositionList[i].subList(IndexMap[i].get(index2), IndexMap[i].get(temp2)));
            else
                second_pos.addAll(PositionList[i].subList(IndexMap[i].get(index2), PositionList[i].size()));
            // build the candidate table
            Candidate_table CT = new Candidate_table();
            for (Integer pos : first_pos) {
                for (int j = mingap; j <= maxgap; j++) {
                    int cand = pos + j + 1;
                    if (cand >= PositionList[i].size())
                        break;
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
        // if support >= minsup, add the candidate to frequent pattern set and store the HBTable in All_HBTS
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

    static void Mining_patterns(Sequence[] SDB, String candidate, List<String>[] frequent) {
        // Mining more_legnth-frequent patterns based on candidate patterns based on 2-frequent patterns HBTables
        int length = candidate.length();
        int support = 0, Length;
        String[] substr = new String[length-1];
        for (int i = 0; i < length-1; i++) {
            substr[i] = candidate.substring(i, i+2);
        }
        // 挖掘3-频繁模式
        if (length == 3) {
            support = Calculate_support(All_HBTS.get(substr[0]), All_HBTS.get(substr[1]), 3);
            if (support >= minsup)
                frequent[length - 1].add(candidate);
        }
        // 挖掘4-频繁模式
        if (length == 4) {
            for (int i = 0; i < SDB_num; i++) {
                Length = SDB[i].length;
                // 先将前两个2-模式拼接成3-模式
                Merge_HBTable(All_HBTS.get(substr[0]).HBTS[i], All_HBTS.get(substr[1]).HBTS[i], temp_HBT);
                // 在最后一次拼接的过程中，边寻找位置组合边完成非重叠位置的筛选
                outerLoop:
                for (int pos : All_HBTS.get(substr[0]).HBTS[i].HBT.keySet()) {
                    for (List<Integer> left_occ : temp_HBT.get(pos)) {
                        if (!All_HBTS.get(substr[length-2]).HBTS[i].HBT.containsKey(left_occ.get(left_occ.size()-1)))
                            continue;
                        for (List<Integer> right_occ : All_HBTS.get(substr[length-2]).HBTS[i].HBT.get(left_occ.get(left_occ.size()-1))) {
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
                for (int j = 0; j < Length; j++)
                    temp_HBT.get(j).clear();
                for (int j = 0; j < length; j++)
                    record[j].clear();
                // 若支持度大于阈值，则提前终止
                if (support >= minsup) {
                    frequent[length - 1].add(candidate);
                    break;
                }
            }
        }
        // 针对5-频繁模式或更长的模式，进行更复杂的挖掘
        if (length >= 5) {
            for (int i = 0; i < SDB_num; i++) {
                Length = SDB[i].length;
                // 先将前两个2-模式拼接成3-模式
                Merge_HBTable(All_HBTS.get(substr[0]).HBTS[i], All_HBTS.get(substr[1]).HBTS[i], temp_HBT);
                // 基于当前的3-模式，进行后续拼接
                for (int l = 2; l < length - 2; l++) {
                    for (int pos : All_HBTS.get(substr[0]).HBTS[i].HBT.keySet()) {
                        List<List<Integer>> temp_list = new ArrayList<>();
                        for (List<Integer> occ : temp_HBT.get(pos)) {
                            List<Integer> temp = new ArrayList<>(occ);
                            temp_list.add(temp);
                        }
                        temp_HBT.get(pos).clear();
                        for (List<Integer> left_occ : temp_list) {
                            if (!All_HBTS.get(substr[l]).HBTS[i].HBT.containsKey(left_occ.get(left_occ.size()-1)))
                                continue;
                            for (List<Integer> right_occ : All_HBTS.get(substr[l]).HBTS[i].HBT.get(left_occ.get(left_occ.size()-1))) {
                                List<Integer> occ = new ArrayList<>(left_occ);
                                occ.add(right_occ.get(0));
                                temp_HBT.get(pos).add(occ);
                            }
                        }
                    }
                }
                // 在最后一次拼接的过程中，边寻找位置组合边完成非重叠位置的筛选
                outerLoop:
                for (int pos : All_HBTS.get(substr[0]).HBTS[i].HBT.keySet()) {
                    for (List<Integer> left_occ : temp_HBT.get(pos)) {
                        if (!All_HBTS.get(substr[length-2]).HBTS[i].HBT.containsKey(left_occ.get(left_occ.size()-1)))
                            continue;
                        for (List<Integer> right_occ : All_HBTS.get(substr[length-2]).HBTS[i].HBT.get(left_occ.get(left_occ.size()-1))) {
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
                for (int j = 0; j < Length; j++)
                    temp_HBT.get(j).clear();
                for (int j = 0; j < length; j++)
                    record[j].clear();
                // 若支持度大于阈值，则提前终止
                if (support >= minsup) {
                    frequent[length - 1].add(candidate);
                    break;
                }
            }
        }
    }

    public static boolean isSubsequence(String shortStr, String longStr) {
        int i = 0, j = 0;
        while (i < shortStr.length() && j < longStr.length()) {
            if (shortStr.charAt(i) == longStr.charAt(j)) {
                i++;
            }
            j++;
        }
        return i == shortStr.length();
    }


    public static void main(String[] args) {

        String[] FileDir = new String[6];
        FileDir[0] = "DataSet/SDB3.txt";
        FileDir[1] = "DataSet/SDB4.txt";
        FileDir[2] = "DataSet/SDB3.txt";
        FileDir[3] = "DataSet/SDB4.txt";

        System.out.println("Welcome to the MNSPM-max!");
        Scanner scanner = new Scanner(System.in);
        int number = scanner.nextInt();
        mingap = scanner.nextInt();
        maxgap = scanner.nextInt();
        double reminsup = scanner.nextDouble();
        scanner.close();
//        int number = Integer.parseInt(args[0]);
//        mingap = Integer.parseInt(args[1]);
//        maxgap = Integer.parseInt(args[2]);
//        double reminsup = Double.parseDouble(args[3]);

        ReadFile(SDB, FileDir[number]);
        minsup = (int) (reminsup * SDB_num);

        List<Integer>[] PositionList = new ArrayList[SDB_num];
        for (int i = 0; i < SDB_num; i++)
            PositionList[i] = new ArrayList<>();
        Map<Integer,Integer>[] IndexList = new HashMap[SDB_num];
        for (int i = 0; i < SDB_num; i++)
            IndexList[i] = new HashMap<>();
        List<String>[] frequent = new ArrayList[maxlen];
        for (int i = 0; i < maxlen; i++)
            frequent[i] = new ArrayList<>();
        for (int i = 0; i < Maxseqlen; i++)
            temp_HBT.put(i, new ArrayList<>());
        for (int i = 0; i < maxlen; i++)
            record[i] = new ArrayList<>();

        MemoryUsageMonitor monitor = new MemoryUsageMonitor();
        monitor.start();
        int begintime = (int) System.currentTimeMillis();

        Build_PositionMap(SDB, PositionList, IndexList);
        Mining_frequent_item(PositionList, IndexList, frequent);
        Prun_PositionMap(PositionList, IndexList, frequent);

        int candidate_num = 0;
        Generate_candidate(candidate, 1, frequent);
        candidate_num += candidate.size();
        for (String pattern : candidate) {
            GapMatch(PositionList, IndexList, pattern, frequent);
        }
        candidate.clear();

        int length = 2;
        while (length < maxlen && !frequent[length - 1].isEmpty()) {
            Generate_candidate(candidate, length, frequent);
            candidate_num += candidate.size();
            // System.out.println(length+1 + "-Candidate Patterns: " + candidate.size());
            length++;
            for (String pattern : candidate) {
                Mining_patterns(SDB, pattern, frequent);
            }
            candidate.clear();
        }

        List<String> maximalPatterns = new ArrayList<>();
        for (int i = 0; i < maxlen; i++) {
            for (String pattern : frequent[i]) {
                boolean isSub = false;
                // 遍历更长的模式，看是否有它的超序列
                for (int j = i + 1; j < maxlen; j++) {
                    for (String superPattern : frequent[j]) {
                        if (isSubsequence(pattern, superPattern)) {
                            isSub = true;
                            break;
                        }
                    }
                    if (isSub) break;
                }
                if (!isSub) {
                    maximalPatterns.add(pattern);
                }
            }
        }

        int endtime = (int) System.currentTimeMillis();
        double maxMemoryUsage = monitor.getMaxMemoryUsage();
        monitor.stopMonitoring();

        System.out.println(FileDir[number]);
        int pattern_num = 0;
        for (int l = 0; l < maxlen; l++) {
            if (frequent[l].isEmpty())
                break;
//            System.out.println(l+1 + "-Frequent Patterns: " + frequent[l].size());
            pattern_num += frequent[l].size();
        }

        System.out.println("The number of candidate patterns: " + candidate_num);
        System.out.println("The number of frequent patterns: " + pattern_num);
        System.out.println("The number of maximal patterns: " + maximalPatterns.size());
        System.out.println("The number of infrequent patterns:" + (candidate_num - pattern_num));

        DecimalFormat df = new DecimalFormat("#.##");
        System.out.println("Time used: " + df.format((endtime - begintime) / 1000.0) + "s");
        System.out.println("Memory used: " + df.format(maxMemoryUsage)+ "MB");

//        List<String> filepaths = new ArrayList<>();
//        filepaths.add("result/gap/SDB1/MNSPM-max.txt");
//        filepaths.add("result/gap/SDB2/MNSPM-max.txt");
//        filepaths.add("result/minsup/SDB1/MNSPM-max.txt");
//        filepaths.add("result/minsup/SDB2/MNSPM-max.txt");
//        try (BufferedWriter writer = new BufferedWriter(new FileWriter(filepaths.get(number), true))) {
//            writer.write(FileDir[number] + "\n");
//            writer.write("The number of candidate patterns: " + candidate_num + "\n");
//            writer.write("The number of frequent patterns: " + pattern_num + "\n");
//            writer.write("Time used: " + df.format((endtime - begintime) / 1000.0) + "s" + "\n");
//            writer.write("Memory used: " + df.format(maxMemoryUsage) + "MB" + "\n");
//        } catch (IOException e) {
//            System.err.println("写入文件时出错: " + e.getMessage());
//        }
    }
}
