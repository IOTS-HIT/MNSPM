// monitor the memory usage of each algorithm
public class MemoryUsageMonitor extends Thread {
    private volatile boolean stopped = false;
    private double maxMemory = 0;

    public void run() {
        while (!stopped) {
            double usedMemory = (Runtime.getRuntime().totalMemory() -  Runtime.getRuntime().freeMemory())
                    / 1024d / 1024d;
            if (usedMemory > maxMemory) {
                maxMemory = usedMemory;
            }
            try {
                Thread.sleep(10); // 每隔1毫秒检测一次
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }

    public double getMaxMemoryUsage() {
        return maxMemory;
    }

    public void stopMonitoring(){
        stopped = true;
    }
}

//import java.lang.management.ManagementFactory;
//import java.lang.management.MemoryMXBean;
//import java.lang.management.MemoryUsage;
//
//public class MemoryUsageMonitor extends Thread {
//    private volatile boolean stopped = false;
//    private double maxMemory = 0;
//    private final MemoryMXBean memoryMXBean = ManagementFactory.getMemoryMXBean();
//
//    public void run() {
//        while (!stopped) {
//            MemoryUsage heapMemoryUsage = memoryMXBean.getHeapMemoryUsage();
//            double usedMemory = heapMemoryUsage.getUsed() / 1024d / 1024d; // 转换为MB
//            if (usedMemory > maxMemory) {
//                maxMemory = usedMemory;
//            }
//            try {
//                Thread.sleep(10); // 每隔10毫秒检测一次
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
//        }
//    }
//
//    public double getMaxMemoryUsage() {
//        return maxMemory;
//    }
//
//    public void stopMonitoring() {
//        stopped = true;
//    }
//}
