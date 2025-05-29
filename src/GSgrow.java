import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

public class GSgrow {

    static int K = 80000;
    static int N = 20000;

    static class SeqDB {
        int id;
        int[] S;
        
        SeqDB(int id, int[] S) {
            this.id = id;
            this.S = S;
        }
    }

    static class Pos {
        List<Integer> pos = new ArrayList<>();
    }

    static class ISupSet {
        int id;
        List<Pos> poset = new ArrayList<>();
    }

    static class FreqPtn {
        int length;
        List<Integer> ptn = new ArrayList<>();
    }

    static class AllFre {
        int length;
        List<FreqPtn> Fre = new ArrayList<>();
    }

    static class SubPtnStruct {
        int start, end;
        int min, max;
    }

    static SeqDB[] sDB = new SeqDB[K];
    static int nn = 0, compnum = 0;
    static List<Integer> eSet = new ArrayList<>();
    static int minsup, NumbS;
    static int mingap, maxgap;
    static int candidate_num = 0;
    static List<FreqPtn> Fre = new ArrayList<>();
    static List<SubPtnStruct> subPtn = new ArrayList<>();
    static List<Integer> S = new ArrayList<>();

    static void readFile(String filename) {
        try (Scanner fileScanner = new Scanner(new File(filename), "UTF-8")) {
            //System.out.println("File opened successfully.");
            int i = 0;
            while (fileScanner.hasNextLine()) {
                String line = fileScanner.nextLine();
                String[] S = line.split("\\s+");
                int[] intArray = java.util.Arrays.stream(S)
                        .filter(s -> !"-1".equals(s) && !"-2".equals(s) && !s.isEmpty())
                        .mapToInt(Integer::parseInt)
                        .toArray();
                sDB[i] = new SeqDB(i, intArray);
//                sDB[i].S = line.toCharArray();
//                sDB[i].id = i + 1;
                i++;
            }
            NumbS = i;
//            System.out.println("The number of sequences is " + NumbS);
        } catch (FileNotFoundException e) {
            System.out.println("Failed to open the file.");
            System.exit(0);
        }
    }

    static void minFreItem() {
        Map<Integer, Integer> counter = new HashMap<>();
        for (int t = 0; t < NumbS; t++) {
            //System.arraycopy(sDB[t].S, 0, S, 0, sDB[t].S.length);
            for (int c : sDB[t].S) {
                S.add(c);
            }
            for (int c : S) {
                counter.put(c, counter.getOrDefault(c, 0) + 1);
            }
            S.clear();
        }
        nn = 0;
        for (Map.Entry<Integer, Integer> entry : counter.entrySet()) {
            if (entry.getValue() >= minsup) {
                eSet.add(entry.getKey());
                nn++;
            }
        }
    }

    static int next(SeqDB seq, int p, int max, int a, int b) {
        for (int i = Math.max(max + 1, a); i <= b; i++) {
            if (seq.S[i] == p) {
                return i;
            }
        }
        return -1;
    }

    static int sizeOneSup(SeqDB[] sDB, List<Integer> P, List<ISupSet> I) {
        int support = 0;
        if (P.size() == 1) {
            for (int i = 0; i < NumbS; i++) {
                for (int l = 0; l < sDB[i].S.length; l++) {
                    if (sDB[i].S[l] == P.get(P.size() - 1)) {
                        while (i >= I.size()) {
                            I.add(new ISupSet());
                        }
                        ISupSet isupSet = I.get(i);
                        isupSet.id = sDB[i].id;
                        Pos pos = new Pos();
                        pos.pos.add(l);
                        isupSet.poset.add(pos);
                        support++;
                    }
                }
            }
        }
        return support;
    }

    static int insgrowGap(SeqDB[] sDB, List<SubPtnStruct> subPtn, List<ISupSet> I) {
        int support = 0;
        int p = subPtn.get(subPtn.size() - 1).end;
        //IPLUS.addAll(I);
        int Size = I.size();
        for (int i = 0; i < Size; i++) {
            int size = I.get(i).poset.size();
            outerloop:
            for (int j = 0; j < size; j++) {
                //System.out.println("i = " + i + " j = " + j);
                Pos apos = I.get(i).poset.get(j);
                int len = apos.pos.size();
                if (len == subPtn.size()) {
                    int flag = 1;
                    for (int k = 0; k < len; k++) {
                        if (sDB[i].S[apos.pos.get(k)] != subPtn.get(k).start) {
                            flag = 0;
                            break;
                        }
                    }
                    if (flag == 0)
                        continue;

                    int max = apos.pos.get(len - 1);
                    int a = max + subPtn.get(subPtn.size() - 1).min + 1;
                    int b = Math.min(max + subPtn.get(subPtn.size() - 1).max + 1, sDB[i].S.length - 1);
                    if (a > sDB[i].S.length) continue;
                    if (max <= b) {
                        int l = next(sDB[i], p, max, a, b);
                        while (l != -1 && l <= b) {
                            for (int m = I.get(i).poset.size()-1; m >= size-1; m--) {
                                if (I.get(i).poset.get(m).pos.size() > len && I.get(i).poset.get(m).pos.get(len) == l) { // 该位置已经被前面的occ用过
                                    l = next(sDB[i], p, l, a, b);
                                    break;
                                }
                                if (m == size-1) {
                                    Pos newPos = new Pos();
                                    newPos.pos.addAll(apos.pos);
                                    newPos.pos.add(l);
                                    I.get(i).poset.add(newPos);
                                    support++;
                                    continue outerloop;
                                }
                            }
                        }
                    }
                }
            }
            for (int k = size-1; k >= 0; k--) {
                if (I.get(i).poset.get(k).pos.size() > subPtn.size() + 1) {
                    I.get(i).poset.remove(k);
                }
            }
        }
        return support;
    }

    static void mineFre(int support, SeqDB[] sDB, List<Integer> P, List<ISupSet> I) {
         candidate_num++;
        if (support >= minsup) {
            //System.out.print(" " + P.toString());
            Fre.add(new FreqPtn());
            FreqPtn fp = Fre.get(Fre.size() - 1);
            fp.length = P.size();
            fp.ptn = new ArrayList<>(P);
            for (int e : eSet) {
                List<Integer> Q = new ArrayList<>(P);
                Q.add(e);
                subPtn.clear();
                for (int k = 0; k < Q.size() - 1; k++) {
                    SubPtnStruct sps = new SubPtnStruct();
                    sps.start = Q.get(k);
                    sps.min = mingap;
                    sps.max = maxgap;
                    sps.end = Q.get(k + 1);
                    subPtn.add(sps);
                    //System.out.print(Q.get(k));
                }
                //System.out.print(Q.get(Q.size() - 1) + " ");
                compnum++;
                //List<ISupSet> IPLUS = new ArrayList<>();
                int newsupport = insgrowGap(sDB, subPtn, I);
                mineFre(newsupport, sDB, Q, I);
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

        System.out.println("Welcome to GSgrow!");
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

        readFile(FileDir[number]);
        minsup = (int) (reminsup * NumbS);
//        System.out.println("The minimum support is: " + minsup);

        MemoryUsageMonitor monitor = new MemoryUsageMonitor();
        monitor.start();
        long starttime = System.currentTimeMillis();

        minFreItem();
//        System.out.println("The number of frequent items: " + nn);
        Fre.clear();

        for (int e = 0; e < nn; e++) {
            List<Integer> P = new ArrayList<>();
            P.add(eSet.get(e));
            List<ISupSet> I = new ArrayList<>();
            int support = sizeOneSup(sDB, P, I);
            mineFre(support, sDB, P, I);
        }

        double maxMemoryUsage = monitor.getMaxMemoryUsage();
        monitor.stopMonitoring();
        long endtime = System.currentTimeMillis();

        System.out.println(FileDir[number]);
        System.out.println("The number of candidate patterns: " + candidate_num);
        System.out.println("The number of frequent patterns: " + Fre.size());
        System.out.println("The number of infrequent patterns: " + (candidate_num - Fre.size()));

        DecimalFormat df = new DecimalFormat("#.##");
        System.out.println("Time used: " + df.format((endtime - starttime)/1000.0) + "s");
        System.out.println("Memory used: " + df.format(maxMemoryUsage) + "MB");

    }
}
