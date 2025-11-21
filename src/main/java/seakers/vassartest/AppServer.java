package seakers.vassartest;

import com.google.gson.*;
import com.google.gson.reflect.TypeToken;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;
import seakers.vassar.evaluation.ArchitectureEvaluationManager;
import seakers.vassar.problems.Assigning.ArchitectureEvaluator;
import seakers.vassar.problems.Assigning.GigaAssigningParams;
import seakers.vassartest.search.problems.Assigning.AssigningProblem;
import seakers.vassartest.search.problems.Assigning.GigaArchitecture;

import java.io.*;
import java.lang.reflect.Type;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantReadWriteLock;

public class AppServer {

    private static volatile boolean INITIALIZED = false;
    private static GigaAssigningParams params;
    private static ArchitectureEvaluationManager evaluationManager;
    private static AssigningProblem problem;
    private static ExecutorService evalExecutor;

    // Guards read access during evaluation vs write access during (re)initialization
    private static final ReentrantReadWriteLock RW = new ReentrantReadWriteLock();

    private static final Gson gson = new GsonBuilder().serializeNulls().create();

    // ---- cache (bitString -> evaluated result) ----
    private static final int MAX_CACHE_SIZE =
            getEnvInt("EVAL_CACHE_SIZE", 10_000); // configurable upper bound

    // Each entry stores the science, cost, and the elapsedMs observed during the
    // original evaluation, plus a timestamp (handy for troubleshooting or TTLs later)
    private static final class CacheEntry {
        final double science;
        final double cost;
        final long elapsedMs;
        final long createdEpochMs;
        CacheEntry(double science, double cost, long elapsedMs) {
            this.science = science;
            this.cost = cost;
            this.elapsedMs = elapsedMs;
            this.createdEpochMs = System.currentTimeMillis();
        }
    }

    // A simple LRU implemented with LinkedHashMap(accessOrder = true)
    // Wrapped in a synchronized Map for thread-safety.
    private static final Map<String, CacheEntry> CACHE =
            Collections.synchronizedMap(new LinkedHashMap<String, CacheEntry>(256, 0.75f, true) {
                @Override
                protected boolean removeEldestEntry(Map.Entry<String, CacheEntry> eldest) {
                    return size() > MAX_CACHE_SIZE;
                }
            });

    // Optional: deduplicate concurrent evaluations of the same bitstring
    // so we don't schedule duplicates while the first is in-flight.
    private static final ConcurrentHashMap<String, Future<EvalResult>> INFLIGHT = new ConcurrentHashMap<>();

    public static void main(String[] args) throws Exception {
        initOnce();

        int port = getEnvInt("PORT", 8080);
        HttpServer server = HttpServer.create(new InetSocketAddress(port), 0);
        server.createContext("/healthz", exchange -> respondJson(exchange, 200, Collections.singletonMap("status","ok")));
        server.createContext("/.well-known/ready", exchange -> respondPlain(exchange, 200, "ready"));
        server.createContext("/evaluate", new EvaluateHandler());
        server.createContext("/initialize", new InitializeHandler()); // <-- NEW
        // Optional: simple cache stats endpoint
        server.createContext("/cachez", exchange -> {
            Map<String, Object> stats = new LinkedHashMap<>();
            synchronized (CACHE) {
                stats.put("size", CACHE.size());
                stats.put("maxSize", MAX_CACHE_SIZE);
            }
            respondJson(exchange, 200, stats);
        });
        server.setExecutor(Executors.newCachedThreadPool());
        server.start();
        System.out.println("Evaluation server listening on port " + port);
    }

    private static synchronized void initOnce() {
        if (INITIALIZED) return;

        int orekitThreads = getEnvInt("OREKIT_THREADS", 1);
        String resourcesPath = getEnvStr("RESOURCES_PATH", "/app/VASSAR_resources");
        int numCpus = getEnvInt("EVAL_CPUS", 1);

        params = new GigaAssigningParams(resourcesPath, "FUZZY-CASES", "test", "normal", orekitThreads);
        ArchitectureEvaluator evaluator = new ArchitectureEvaluator();
        evaluationManager = new ArchitectureEvaluationManager(params, evaluator);
        evaluationManager.init(numCpus);

        problem = new AssigningProblem(new int[]{1}, "GigaProblem", evaluationManager, params);

        evalExecutor = Executors.newSingleThreadExecutor(new ThreadFactory() {
            @Override public Thread newThread(Runnable r) {
                Thread t = new Thread(r, "eval-exec");
                t.setDaemon(true);
                return t;
            }
        });

        INITIALIZED = true;
    }

    // --- made static and write-locked ---
    private static void initPanelWeights(HashMap<String, Double> panelWeights) {
        RW.writeLock().lock();
        try {
            int orekitThreads = getEnvInt("OREKIT_THREADS", 1);
            String resourcesPath = getEnvStr("RESOURCES_PATH", "/app/VASSAR_resources");
            int numCpus = getEnvInt("EVAL_CPUS", 1);

            GigaAssigningParams newParams = new GigaAssigningParams(resourcesPath, "FUZZY-CASES", "test", "normal", orekitThreads);
            if (panelWeights != null) {
                newParams.setPanelWeightMap(panelWeights);
            }

            ArchitectureEvaluator evaluator = new ArchitectureEvaluator();
            ArchitectureEvaluationManager newManager = new ArchitectureEvaluationManager(newParams, evaluator);
            newManager.init(numCpus);

            // Atomically swap shared state
            params = newParams;
            evaluationManager = newManager;
            problem = new AssigningProblem(new int[]{1}, "GigaProblem", evaluationManager, params);

            // Invalidate caches because the scoring surface changed
            CACHE.clear();
            INFLIGHT.clear();
        } finally {
            RW.writeLock().unlock();
        }
    }

    static class EvaluateHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange exchange) throws IOException {
            if (!"POST".equalsIgnoreCase(exchange.getRequestMethod())) {
                respondJson(exchange, 405, Collections.singletonMap("error", "Use POST"));
                return;
            }
            if (!contentTypeIsJson(exchange)) {
                respondJson(exchange, 415, Collections.singletonMap("error", "Content-Type must be application/json"));
                return;
            }

            String body = new String(exchange.getRequestBody().readAllBytes(), StandardCharsets.UTF_8);

            try {
                Type mapType = new TypeToken<Map<String, List<String>>>(){}.getType();
                JsonObject root = JsonParser.parseString(body).getAsJsonObject();

                // Accept either: { "bitString": "..." } OR { "design": {orbit: [instruments]} }
                String bitString = null;
                if (root.has("bitString") && root.get("bitString").isJsonPrimitive()) {
                    bitString = root.get("bitString").getAsString();
                } else {
                    if (!root.has("design") || !root.get("design").isJsonObject()) {
                        respondJson(exchange, 400, Collections.singletonMap("error", "Provide either 'bitString' or a 'design' object (orbit -> [instruments])"));
                        return;
                    }
                    Map<String, List<String>> designMap = gson.fromJson(root.get("design"), mapType);

                    HashMap<String, ArrayList<String>> hm = new HashMap<>();
                    for (Map.Entry<String, List<String>> e : designMap.entrySet()) {
                        hm.put(e.getKey(), new ArrayList<String>(e.getValue()));
                    }
                    bitString = params.getBitString(hm);
                }

                // Check to see if re-initialization is requested in the evaluation
                if (root.has("panelWeights") && root.get("panelWeights").isJsonObject()) {
                    JsonObject pwObj = root.get("panelWeights").getAsJsonObject();
                    HashMap<String, Double> panelWeights = new HashMap<>();
                    for (Map.Entry<String, JsonElement> e : pwObj.entrySet()) {
                        panelWeights.put(e.getKey(), e.getValue().getAsDouble());
                    }
                    initPanelWeights(panelWeights);
                }

                final String finalBitString = bitString;
                long t0 = System.nanoTime();

                // ---- CACHE USAGE ----
                // 1) Fast path: cache hit
                CacheEntry cached;
                synchronized (CACHE) { // LinkedHashMap LRU needs external synchronization
                    cached = CACHE.get(finalBitString);
                }
                if (cached != null) {
                    Map<String, Object> response = new LinkedHashMap<>();
                    response.put("science", cached.science);
                    response.put("cost", cached.cost);
                    response.put("bitString", finalBitString);
                    response.put("elapsedMs", cached.elapsedMs); // original eval latency
                    response.put("cached", true);
                    respondJson(exchange, 200, response);
                    return;
                }

                // 2) Cache miss: deduplicate concurrent evaluations for the same bitString
                Future<EvalResult> fut = INFLIGHT.computeIfAbsent(finalBitString, bs ->
                        evalExecutor.submit(new Callable<EvalResult>() {
                            @Override public EvalResult call() throws Exception {
                                RW.readLock().lock(); // prevent re-init during evaluation
                                try {
                                    GigaArchitecture arch = new GigaArchitecture(bs);
                                    problem.evaluate(arch);
                                    double science = -1.0 * arch.getObjective(0);
                                    double cost = arch.getObjective(1);
                                    return new EvalResult(bs, science, cost);
                                } finally {
                                    RW.readLock().unlock();
                                }
                            }
                        })
                );

                // Wait for the (possibly shared) result
                EvalResult result;
                try {
                    result = fut.get();
                } finally {
                    // Ensure INFLIGHT is cleaned so future identical requests can re-enter
                    INFLIGHT.remove(finalBitString, fut);
                }

                long elapsedMs = TimeUnit.NANOSECONDS.toMillis(System.nanoTime() - t0);

                // 3) Store in cache for future requests using the original elapsed
                CacheEntry toCache = new CacheEntry(result.science, result.cost, elapsedMs);
                synchronized (CACHE) {
                    CACHE.put(finalBitString, toCache);
                }

                Map<String, Object> response = new LinkedHashMap<String, Object>();
                response.put("science", result.science);
                response.put("cost", result.cost);
                response.put("bitString", result.bitString);
                response.put("elapsedMs", elapsedMs);
                response.put("cached", false);

                respondJson(exchange, 200, response);

            } catch (RejectedExecutionException rex) {
                Map<String, Object> resp = new LinkedHashMap<String, Object>();
                resp.put("error", "Server busy, retry");
                resp.put("details", rex.getMessage());
                respondJson(exchange, 503, resp);
            } catch (Exception ex) {
                ex.printStackTrace();
                Map<String, Object> resp = new LinkedHashMap<String, Object>();
                resp.put("error", "Evaluation failed");
                resp.put("details", ex.getMessage());
                respondJson(exchange, 500, resp);
            }
        }
    }

    // --- NEW: /initialize handler ---
    static class InitializeHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange exchange) throws IOException {
            if (!"POST".equalsIgnoreCase(exchange.getRequestMethod())) {
                respondJson(exchange, 405, Collections.singletonMap("error", "Use POST"));
                return;
            }
            if (!contentTypeIsJson(exchange)) {
                respondJson(exchange, 415, Collections.singletonMap("error", "Content-Type must be application/json"));
                return;
            }

            String body = new String(exchange.getRequestBody().readAllBytes(), StandardCharsets.UTF_8);
            try {
                JsonObject root = JsonParser.parseString(body).getAsJsonObject();
                if (!root.has("panelWeights") || !root.get("panelWeights").isJsonObject()) {
                    respondJson(exchange, 400, Collections.singletonMap("error", "Provide 'panelWeights' object"));
                    return;
                }

                JsonObject pwObj = root.getAsJsonObject("panelWeights");
                HashMap<String, Double> panelWeights = new HashMap<>();
                for (Map.Entry<String, JsonElement> e : pwObj.entrySet()) {
                    if (!e.getValue().isJsonPrimitive() || !e.getValue().getAsJsonPrimitive().isNumber()) {
                        respondJson(exchange, 400, Collections.singletonMap("error", "All panelWeights values must be numbers"));
                        return;
                    }
                    panelWeights.put(e.getKey(), e.getValue().getAsDouble());
                }

                initPanelWeights(panelWeights);

                Map<String, Object> resp = new LinkedHashMap<>();
                resp.put("status", "ok");
                resp.put("message", "Initialization completed with provided panel weights.");
                resp.put("numWeights", panelWeights.size());
                respondJson(exchange, 200, resp);

            } catch (Exception ex) {
                ex.printStackTrace();
                Map<String, Object> resp = new LinkedHashMap<>();
                resp.put("error", "Initialization failed");
                resp.put("details", ex.getMessage());
                respondJson(exchange, 500, resp);
            }
        }
    }

    // ---- helpers ----
    private static boolean contentTypeIsJson(HttpExchange exchange) {
        List<String> ct = exchange.getRequestHeaders().get("Content-Type");
        if (ct == null || ct.isEmpty()) return false;
        for (String s : ct) {
            if (s != null && s.toLowerCase().contains("application/json")) return true;
        }
        return false;
    }

    private static void respondJson(HttpExchange exchange, int status, Object obj) throws IOException {
        byte[] payload = gson.toJson(obj).getBytes(StandardCharsets.UTF_8);
        exchange.getResponseHeaders().add("Content-Type", "application/json; charset=utf-8");
        exchange.getResponseHeaders().add("Access-Control-Allow-Origin", "*");
        exchange.sendResponseHeaders(status, payload.length);
        OutputStream os = exchange.getResponseBody();
        try { os.write(payload); }
        finally { os.close(); }
    }

    private static void respondPlain(HttpExchange exchange, int status, String text) throws IOException {
        byte[] payload = text.getBytes(StandardCharsets.UTF_8);
        exchange.getResponseHeaders().add("Content-Type", "text/plain; charset=utf-8");
        exchange.sendResponseHeaders(status, payload.length);
        OutputStream os = exchange.getResponseBody();
        try { os.write(payload); }
        finally { os.close(); }
    }

    private static String getEnvStr(String key, String fallback) {
        String v = System.getenv(key);
        if (v == null || v.trim().isEmpty()) return fallback;
        return v;
    }

    private static int getEnvInt(String key, int fallback) {
        try {
            String v = System.getenv(key);
            return (v == null || v.trim().isEmpty()) ? fallback : Integer.parseInt(v.trim());
        } catch (Exception e) {
            return fallback;
        }
    }

    // simple POJO instead of 'record'
    private static class EvalResult {
        final String bitString;
        final double science;
        final double cost;
        EvalResult(String bitString, double science, double cost) {
            this.bitString = bitString;
            this.science = science;
            this.cost = cost;
        }
    }
}
