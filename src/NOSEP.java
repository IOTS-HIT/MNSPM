import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class NOSEP {

    static int K = 80000; // The sequence number of sequence database
    static int M = 100; // The length of pattern
    static int N = 1000; // The length of sequence

    static class Node {
        //public boolean toleave;
        int name; // The corresponding position of node in sequence
        int minLeave, maxLeave; // The position of maximum and minimum leave nodes
        List<Integer> parent; // The position set of parents
        List<Integer> children; // The position set of children
        boolean used; // true if has used, false if has not used
        boolean toLeave; // true if can reach leaves, false if not

        Node() {
            parent = new ArrayList<>();
            children = new ArrayList<>();
        }
    }
    static class SeqDB {
        int id;
        int[] S;

         SeqDB(int id, int[] S) {
             this.id = id;
             this.S = S;
         }
    }

    static class Pattern {
        int[] pattern;
        int length;

        Pattern(int[] pattern, int length) {
            this.pattern = pattern;
            this.length = length;
        }
    }

    static class Occurrence {
        List<Integer> position;

        Occurrence() {
            position = new ArrayList<>();
        }
    }

    static class SubPattern {
        int start, end;
        int min, max;
    }

    static int store, ptnLen;
    static int mingap, maxgap; // gap constraint
    static int[] S = new int[N]; // sequence
    static int NumbS, minsup;
    static SubPattern[] subPtn = new SubPattern[M]; // pattern p[i]
    static List<Pattern>[] freArr = new ArrayList[M]; // store frequent patterns
    static List<Pattern> candidate = new ArrayList<>(); // store candidate patterns

    static void readFile(SeqDB[] sDB, String filename) {
        try (Scanner fileScanner = new Scanner(new File(filename), "UTF-8")) {
//            System.out.println("File opened successfully.");
            int i = 0;
            while (fileScanner.hasNextLine() && i < K) {
                String line = fileScanner.nextLine();
                String[] S = line.split("\\s+");
                int[] intArray = java.util.Arrays.stream(S)
                        .filter(s -> !"-1".equals(s) && !"-2".equals(s) && !s.isEmpty())
                        .mapToInt(Integer::parseInt)
                        .toArray();
                sDB[i] = new SeqDB(i, intArray);
                i++;
            }
            NumbS = i;
        } catch (FileNotFoundException e) {
            System.out.println("Failed to open the file.");
            System.exit(0);
        }
    }

    static void minFreItem(SeqDB[] sDB) {
        // Initialize freArr
        for (int i = 0; i < M; i++) {
            freArr[i] = new ArrayList<>();
        }
        Map<Integer, Integer> counter = new HashMap<>();
        String mine;
        for (int t = 0; t < NumbS; t++) {
            int[] sequence = sDB[t].S;
            for (int c : sequence) {
                counter.put(c, counter.getOrDefault(c, 0) + 1);
            }
        }
        Iterator<Map.Entry<Integer, Integer>> iterator = counter.entrySet().iterator();
        while (iterator.hasNext()) {
            Map.Entry<Integer, Integer> entry = iterator.next();
            if (entry.getValue() < minsup) {
                iterator.remove();
            } else {
                Pattern p = new Pattern(new int[]{entry.getKey()}, 1);
                freArr[0].add(p);
            }
        }
    }

    static void genCandidate(int level) {
        if (level == 1) {
            for (Pattern p : freArr[0]) {
                for (Pattern q : freArr[0]) {
                    Pattern super_pattern = new Pattern(new int[]{p.pattern[0], q.pattern[0]}, 2);
                    candidate.add(super_pattern);
                }
            }
        }
        else {
            for (Pattern p : freArr[level - 1]) {
                for (Pattern q : freArr[level - 1]) {
                    int flag = 0;
                    for (int i = 0; i < level - 1; i++){
                        if (p.pattern[i] != q.pattern[i+1]){
                            flag = 1;
                            break;
                        }
                    }
                    if (flag == 0){
                        int[] pattern = Arrays.copyOf(q.pattern, level+1);
                        pattern[level] = p.pattern[level-1];
                        Pattern super_pattern = new Pattern(pattern, level+1);
                        candidate.add(super_pattern);
                    }
                }
            }
        }
    }

    static void dealRange(Pattern pattern) {
        ptnLen = 0;
        //Arrays.fill(subPtn, new SubPattern());
        for (int i = 0; i < pattern.length; i++) {
            subPtn[i] = new SubPattern();
        }
        int[] p = pattern.pattern;
        if (p.length == 1) {
            subPtn[ptnLen].start = p[0];
            subPtn[ptnLen].max = subPtn[ptnLen].min = 0;
        }
        for (int i = 0; i < p.length - 1; i++) {
            subPtn[ptnLen].start = p[i];
            subPtn[ptnLen].end = p[i + 1];
            subPtn[ptnLen].min = mingap;
            subPtn[ptnLen].max = maxgap;
            ptnLen++;
        }
    }

    static void createNettree(List<Node>[] nettree) {
        //for (int i = 0; i < ptnLen + 1; i++)
        //nettree[i] = new ArrayList<>();
        int[] start = new int[ptnLen + 1];
        Arrays.fill(start, 0);
        int len = S.length;
        for (int i = 0; i < len; i++) {
            Node anode = new Node();
            anode.name = i;
            anode.parent = new ArrayList<>();
            anode.children = new ArrayList<>();
            anode.maxLeave = anode.name;
            anode.minLeave = anode.name;
            anode.used = false;
            if (subPtn[0].start == S[i]) {
                anode.toLeave = true;
                nettree[0].add(anode);
            }
            for (int j = 0; j < ptnLen; j++) {
                if (subPtn[j].end == S[i]) {
                    int prevLen = nettree[j].size();
                    if (prevLen == 0) break;
                    for (int k = start[j]; k < prevLen; k++) {
                        if (i - nettree[j].get(k).name - 1 > subPtn[j].max) {
                            start[j]++;
                        }
                    }
                    if (i - nettree[j].get(prevLen - 1).name - 1 > subPtn[j].max) continue;
                    if (i - nettree[j].get(start[j]).name - 1 < subPtn[j].min) continue;
                    Node newNode = new Node();
                    newNode.name = i;
                    newNode.parent = new ArrayList<>();
                    newNode.children = new ArrayList<>();
                    newNode.maxLeave = anode.name;
                    newNode.minLeave = anode.name;
                    newNode.used = false;
                    newNode.toLeave = true;
                    nettree[j + 1].add(newNode);
                    for (int k = start[j]; k < prevLen; k++) {
                        if (i - nettree[j].get(k).name - 1 < subPtn[j].min) break;
                        int nc = nettree[j].get(k).children.size();
                        nettree[j].get(k).children.add(nettree[j + 1].size() - 1);
                        nettree[j + 1].get(nettree[j + 1].size() - 1).parent.add(k);
                    }
                }
            }
        }
    }

    static void updateNettree(List<Node>[] nettree) {
        for (int i = ptnLen - 1; i >= 0; i--) {
            for (int j = nettree[i].size() - 1; j >= 0; j--) {
                boolean flag = true;
                int size = nettree[i].get(j).children.size();
                for (int k = 0; k < size; k++) {
                    int child = nettree[i].get(j).children.get(k);
                    if (k == 0) {
                        nettree[i].get(j).minLeave = nettree[i + 1].get(child).minLeave;
                    }
                    if (k == size - 1) {
                        nettree[i].get(j).maxLeave = nettree[i + 1].get(child).maxLeave;
                    }
                    if (!nettree[i + 1].get(child).used) {
                        flag = false;
                    }
                }
                nettree[i].get(j).used = flag;
                if (flag) {
                    nettree[i].get(j).maxLeave = nettree[i].get(j).name;
                    nettree[i].get(j).minLeave = nettree[i].get(j).name;
                    nettree[i].get(j).toLeave = false;
                }
            }
        }
    }

    static void updateNettreePC(List<Node>[] nettree, Occurrence occin) {
        for (int level = ptnLen; level > 0; level--) {
            int position = occin.position.get(level);
            int num = nettree[level].size();
            for (; position < num; position++) {
                if (!nettree[level].get(position).used) break;
                int len = nettree[level].get(position).parent.size();
                for (int i = 0; i < len; i++) {
                    int parent = nettree[level].get(position).parent.get(i);
                    int cs = nettree[level - 1].get(parent).children.size();
                    if (nettree[level - 1].get(parent).used) continue;
                    if (cs == 1) {
                        nettree[level - 1].get(parent).used = true;
                        nettree[level - 1].get(parent).toLeave = false;
                    } else {
                        int kk;
                        for (kk = 0; kk < cs; kk++) {
                            int child = nettree[level - 1].get(parent).children.get(kk);
                            if (!nettree[level].get(child).used) break;
                        }
                        if (kk == cs) {
                            nettree[level - 1].get(parent).used = true;
                            nettree[level - 1].get(parent).toLeave = false;
                        }
                    }
                }
            }
        }
    }

    static void nonOverLength(int rest, List<Node>[] nettree) {
        createNettree(nettree);
        updateNettree(nettree);
        store = 0;
        int num = -1;
        for (Node node : nettree[0]) {
            num++;
            if (!node.toLeave) continue;
            node.used = true;
            node.toLeave = false;
            Occurrence occin = new Occurrence();
            //occin.position.add(0, node.name);
            occin.position.add(0, num);
            node.used = true;
            node.toLeave = false;
            int j;
            for (j = 1; j < ptnLen + 1; j++) {
                int parent = occin.position.get(j - 1);
                //int cs = nettree[j - 1].get(num).children.size();
                int cs = nettree[j - 1].get(parent).children.size();
                int t;
                for (t = 0; t < cs; t++) {
                    //int child = nettree[j - 1].get(num).children.get(t);
                    int child = nettree[j - 1].get(parent).children.get(t);
                    if (!nettree[j].get(child).used) {
                        occin.position.add(j, child);
                        //occin.position.add(j, nettree[j].get(child).name);
                        nettree[j].get(child).used = true;
                        nettree[j].get(child).toLeave = false;
                        break;
                    }
                }
                if (t == cs) {
                    for (int kk = 0; kk < j; kk++) {
                        int pos = occin.position.get(kk);
                        nettree[kk].get(pos).used = false;
                        nettree[kk].get(pos).toLeave = true;
                    }
                    break;
                }
            }
            if (j == ptnLen + 1) {
                store++;
                if (store > rest) break;
                updateNettreePC(nettree, occin);
            }
            occin.position.clear();
        }
    }

    static int netGap(Pattern p, int rest) {
        dealRange(p);
        if (ptnLen + 1 > S.length) return 0;
        List<Node>[] nettree = new ArrayList[ptnLen + 1];;
        for (int i = 0; i < ptnLen + 1; i++) {
            nettree[i] = new ArrayList<>();
        }
        nonOverLength(rest, nettree);
        return store;
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

        System.out.println("Welcome to NOSEP!");
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

        SeqDB[] sDB = new SeqDB[K];
        readFile(sDB, FileDir[number]);
        minsup = (int) (reminsup * NumbS);
//        System.out.println("The minimum support is: " + minsup);

        MemoryUsageMonitor monitor = new MemoryUsageMonitor();
        monitor.start();

        int begintime = (int) System.currentTimeMillis();
        minFreItem(sDB);
//        System.out.println("The number of frequent items: " + freArr[0].size());

        int fLevel = 1;
        int candidate_num = 0;

        genCandidate(fLevel);
//        System.out.println("The number of 2-candidate patterns: " + candidate.size());
        candidate_num += candidate.size();

        while (!candidate.isEmpty()) {
            int num = 0;
            for (Pattern p : candidate) {
                int occnum = 0;
                int rest = 0;
                for (int t = 0; t < NumbS; t++) {
                    rest = minsup - occnum;
                    if (sDB[t].S.length > 0) {
                        S = sDB[t].S;
                        occnum += netGap(p, rest);
                    }
                    if (occnum >= minsup) {
                        freArr[fLevel].add(p);
                        break;
                    }
                }
            }
            fLevel++;
            candidate.clear();
            genCandidate(fLevel);
            candidate_num += candidate.size();
        }

        int endtime = (int) System.currentTimeMillis();

        double maxMemory = monitor.getMaxMemoryUsage();
        monitor.stopMonitoring();

        int frenum = 0;
        for (int i = 0; i < fLevel; i++) {
            if (freArr[i].isEmpty())
                break;
//            System.out.println(i+1 + "-Frequent Patterns: " + freArr[i].size());
            frenum += freArr[i].size();
        }

        System.out.println(FileDir[number]);
        System.out.println("The number of candidate patterns: " + candidate_num);
        System.out.println("The number of frequent patterns: " + frenum);
        System.out.println("The number of infrequent patterns: " + (candidate_num - frenum));

        DecimalFormat df = new DecimalFormat("#.##");
        System.out.println("Time used: " + df.format((endtime - begintime) / 1000.0) + "s");
        System.out.println("Memory used: " + df.format(maxMemory) + "MB");

    }
}
