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