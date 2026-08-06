/**
 * @file generator.c
 * @brief C Program for Polymer Network Generation via Strict Sculpting (Serial Version with Logging).
 *
 * ROLE: this is the standalone searcher. It runs on its own, without
 * Python, and is the tool for long exhaustive searches over many trials.
 * The pure-Python port in topon/topology/generator_python.py is the
 * separate quick path for in-process generation of likely networks.
 * The two are deliberately independent programs, not a library and a
 * wrapper: nothing here is called from Python, and nothing here should
 * grow a Python binding.
 *
 * PROVENANCE: vendored 2026-08-05 from generator_serial_debug11.c
 * (md5 e7631f4bbcb963d50c382721de3b3c18, dated 2025-11-03), the version
 * the shipped generator.exe was built from and the one the Python port
 * mirrors.
 *
 * A later variant exists in the archive (md5 83d7f9d3, 2026-02-27, under
 * experiments/pruning_research/pruning_algorithm_math*) which replaces
 * the per-degree count check in is_move_safe with a cumulative one. It is
 * NOT used here: measured across six standard SC configurations it
 * sculpts 1/6 where this version sculpts 6/6, failing whenever max_func
 * is below the lattice coordination. Treat it as an open experiment.
 *
 * Build:  gcc -O2 -o generator.exe generator.c -lm
 *
 * Lattice construction and the .nodes/.edges format are shared surface
 * with the Python port and must be changed in both; see
 * tests/unit/topology/test_c_generator.py. The sculpting search itself
 * is this program's own business.
 *
 * KNOWN DIVERGENCES from the Python port (pre-existing, not fixed here):
 *   - Periodicity: this file honours p_dims per axis; the Python
 *     builders always wrap.
 *   - The degree<=2 guard in the sculpting stages is gated on
 *     is_sc_lattice here, but applied unconditionally in Python.
 *
 * @details
 * This program simulates the creation of a polymer network with a specific target
 * degree distribution. It employs a **Strict Sculpting Model** to rigorously avoid
 * unintended "collateral damage," especially the creation of nodes with degrees
 * that have a target count of zero.
 *
 * MODIFIED: Now supports SC, BCC, and FCC initial lattice generation.
 * MODIFIED: Now supports an 'e:N' argument to target a specific *total edge count*.
 *
 * @author Ahmet Burak Yildirim
 * @date November 1, 2025
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <math.h> // For sqrt
#include <sys/stat.h> // For mkdir
#ifdef _WIN32
#include <direct.h> // For _mkdir on Windows
#include <process.h> // For _getpid
#define mkdir(dir, mode) _mkdir(dir)
#define topon_getpid _getpid
#else
#include <unistd.h> // For getpid
#define topon_getpid getpid
#endif


// --- Structs and Type Definitions ---

// --- NEW: Coordinate struct ---
typedef struct Coord {
    double x, y, z;
} Coord;

typedef enum NodeStatus {
    ACTIVE,
    IS_DEGREE_0,
    IS_DEGREE_1
} NodeStatus;

typedef enum MoveType {
    SET_D0,
    SET_D1,
    REMOVE_EDGE
} MoveType;

typedef struct MoveLog {
    MoveType type;
    int u;
    int v; // Second node for edges, -1 for node operations
} MoveLog;

typedef struct Edge {
    int u;
    int v;
} Edge;

typedef struct AdjListNode {
    int dest;
    struct AdjListNode* next;
} AdjListNode;

typedef struct AdjList {
    AdjListNode *head;
} AdjList;

typedef struct Graph {
    int V;
    AdjList* array;
    int* degrees;
    Coord* coords; // --- MODIFIED: Added coordinates array ---
} Graph;

typedef struct UnionFind {
    int* parent;
    int n;
} UnionFind;

// --- Forward declarations ---
// --- MODIFIED: Signatures updated for new 'e:N' logic ---
Graph* run_single_trial(Graph* base_graph, int max_func, const int* target_counts, int target_edge_count, long long trial_num, int extensive_logging, const char* dims_str, const char* lattice_type);
void addEdge(Graph* graph, int src, int dest);
void print_distribution(const char* stage_name, Graph* g, const int* target_counts, long long trial_num, long long move_num, int max_func, int extensive_logging);
void save_move_log_to_file(const MoveLog* move_log, long long count, const char* dims_str, long long trial);
int is_move_safe(Graph* g, int u, int v, const int* target_counts, int max_func, int current_stage, long long target_degree_sum, long long current_total_degree_sum);


// --- START: HELPER FUNCTION ---
/**
 * @brief Checks if removing an edge between u and v would violate strict count constraints.
 * @param current_stage The stage calling the function (1=SetD0, 2=SetD1, 3=SetMax, 4=Systematic)
 * @param target_degree_sum The target total degree sum (2 * e_target), or -2 if not set.
 * @param current_total_degree_sum The current total degree sum.
 * @return 1 if the move is safe, 0 if it is forbidden.
 */
// --- MODIFIED: Signature updated for new 'e:N' logic ---
int is_move_safe(Graph* g, int u, int v, const int* target_counts, int max_func, int current_stage, long long target_degree_sum, long long current_total_degree_sum) {
    
    // --- NEW: Target Edge Count Check ---
    // This is the primary "stop" signal when using 'e:N'
    if (current_stage == 4 && target_degree_sum != -1) {
        // If we are at or below the target, no more moves are safe.
        // The degree sum *before* this move is current_total_degree_sum.
        // After this move, it will be current_total_degree_sum - 2.
        if (current_total_degree_sum <= target_degree_sum) {
            return 0; // FORBIDDEN: We have hit or gone below our target edge count.
        }
    }
    // --- END NEW CHECK ---
    
    int u_new_degree = g->degrees[u] - 1;
    int v_new_degree = g->degrees[v] - 1;

    // --- Check 1: Check Neighbor 'v' (The "Victim") ---
    // 'v' is always checked for collateral damage, regardless of stage.

    // 1a. Forbidden Degree (target=0)
    if (v_new_degree >= 0 && target_counts[v_new_degree] == 0) {
        return 0; // Forbidden collateral damage on v
    }

    // 1b. Overshooting
    if (v_new_degree >= 0 && target_counts[v_new_degree] > 0) {
        
        // Case A: d0 or d1. These are "sacred" after Stages 1 & 2.
        // Never overshoot them in *any* stage.
        if (v_new_degree <= 1) { 
            int current_count = 0;
            for (int i = 0; i < g->V; ++i) {
                if (g->degrees[i] == v_new_degree) current_count++;
            }
            if (current_count >= target_counts[v_new_degree]) {
                return 0; // Forbid overshooting d0 or d1
            }
        }
        
        // Case B: d2+. Only block overshooting explicit targets in Stage 4.
        // Stages 1, 2, and 3 are *allowed* to overshoot these.
        else if (current_stage == 4) { 
            int current_count = 0;
            for (int i = 0; i < g->V; ++i) {
                if (g->degrees[i] == v_new_degree) current_count++;
            }
            if (current_count >= target_counts[v_new_degree]) {
                return 0; // Forbid overshooting v's target in Stage 4
            }
        }
        // (If stage 1, 2, or 3, we allow overshooting v for d2+)
    }

    // --- Check 2: Check Actor 'u' ---
    // 'u' is *only* checked for damage in Stage 4.
    // In Stages 1, 2, 3, 'u' is the node we are *trying* to change,
    // so its new state isn't "collateral damage."

    if (current_stage == 4) {
        // 2a. Forbidden Degree (target=0)
        if (u_new_degree >= 0 && target_counts[u_new_degree] == 0) {
            return 0; // Forbid u from becoming a forbidden degree
        }

        // 2b. Overshooting (Only for *explicit* targets)
        if (u_new_degree >= 0 && target_counts[u_new_degree] > 0) {
            int current_count = 0;
            for (int i = 0; i < g->V; ++i) {
                if (g->degrees[i] == u_new_degree) current_count++;
            }
            if (current_count >= target_counts[u_new_degree]) {
                return 0; // Forbid overshooting u's target
            }
        }
    }
    
    // If we passed all checks (e.g., in Stage 1, Check 2 is skipped),
    // the move is safe.
    return 1;
}
// --- END: HELPER FUNCTION ---


// --- Union-Find Data Structure Functions ---

UnionFind* createUnionFind(int n) {
    UnionFind* uf = (UnionFind*)malloc(sizeof(UnionFind));
    uf->parent = (int*)malloc(n * sizeof(int));
    uf->n = n;
    for (int i = 0; i < n; i++) uf->parent[i] = i;
    return uf;
}

int find_set(UnionFind* uf, int i) {
    if (uf->parent[i] == i) return i;
    return uf->parent[i] = find_set(uf, uf->parent[i]);
}

void unite_sets(UnionFind* uf, int a, int b) {
    a = find_set(uf, a);
    b = find_set(uf, b);
    if (a != b) uf->parent[b] = a;
}

void freeUnionFind(UnionFind* uf) {
    if (!uf) return;
    free(uf->parent);
    free(uf);
}

// --- Graph Utility Functions ---

AdjListNode* newAdjListNode(int dest) {
    AdjListNode* newNode = (AdjListNode*)malloc(sizeof(AdjListNode));
    newNode->dest = dest;
    newNode->next = NULL;
    return newNode;
}

// --- MODIFIED: createGraph now allocates space for coordinates ---
Graph* createGraph(int V) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->V = V;
    graph->array = (AdjList*)malloc(V * sizeof(AdjList));
    graph->degrees = (int*)calloc(V, sizeof(int));
    graph->coords = (Coord*)malloc(V * sizeof(Coord)); // Allocate coords
    for (int i = 0; i < V; ++i) {
        graph->array[i].head = NULL;
        graph->coords[i] = (Coord){0.0, 0.0, 0.0}; // Initialize
    }
    return graph;
}

// --- MODIFIED: freeGraph now frees coordinates ---
void freeGraph(Graph* graph) {
    if (!graph) return;
    for (int i = 0; i < graph->V; ++i) {
        AdjListNode* pCrawl = graph->array[i].head;
        while (pCrawl) {
            AdjListNode* temp = pCrawl;
            pCrawl = pCrawl->next;
            free(temp);
        }
    }
    free(graph->array);
    free(graph->degrees);
    free(graph->coords); // Free coords
    free(graph);
}

// --- MODIFIED: copyGraph now copies coordinates ---
Graph* copyGraph(Graph* src_graph) {
    if (!src_graph) return NULL;
    Graph* new_graph = createGraph(src_graph->V);
    // Copy coordinates
    memcpy(new_graph->coords, src_graph->coords, src_graph->V * sizeof(Coord));
    for (int i = 0; i < src_graph->V; i++) {
        AdjListNode* pCrawl = src_graph->array[i].head;
        while(pCrawl) {
            if (i < pCrawl->dest) addEdge(new_graph, i, pCrawl->dest);
            pCrawl = pCrawl->next;
        }
    }
    return new_graph;
}

void addEdge(Graph* graph, int src, int dest) {
    AdjListNode* newNode = newAdjListNode(dest);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
    graph->degrees[src]++;
    newNode = newAdjListNode(src);
    newNode->next = graph->array[dest].head;
    graph->array[dest].head = newNode;
    graph->degrees[dest]++;
}

void removeEdge(Graph* graph, int src, int dest) {
    AdjListNode* pCrawl = graph->array[src].head;
    AdjListNode* prev = NULL;
    while (pCrawl && pCrawl->dest != dest) { prev = pCrawl; pCrawl = pCrawl->next; }
    if (pCrawl) {
        if (prev) prev->next = pCrawl->next; else graph->array[src].head = pCrawl->next;
        free(pCrawl);
        graph->degrees[src]--;
    }
    pCrawl = graph->array[dest].head;
    prev = NULL;
    while (pCrawl && pCrawl->dest != src) { prev = pCrawl; pCrawl = pCrawl->next; }
    if (pCrawl) {
        if (prev) prev->next = pCrawl->next; else graph->array[dest].head = pCrawl->next;
        free(pCrawl);
        graph->degrees[dest]--;
    }
}

// --- Connectivity and File I/O ---

int is_subgraph_connected(Graph* g, const NodeStatus* node_status) {
    if (g->V == 0) return 1;
    UnionFind* uf = createUnionFind(g->V);
    int first_active_node = -1;
    for (int i = 0; i < g->V; i++) {
        if (node_status[i] == ACTIVE) {
            if (first_active_node == -1) first_active_node = i;
            AdjListNode* pCrawl = g->array[i].head;
            while (pCrawl) {
                if (node_status[pCrawl->dest] == ACTIVE) unite_sets(uf, i, pCrawl->dest);
                pCrawl = pCrawl->next;
            }
        }
    }
    if (first_active_node == -1) { freeUnionFind(uf); return 1; }
    int root = find_set(uf, first_active_node);
    int connected = 1;
    for (int i = 0; i < g->V; i++) {
        if (node_status[i] == ACTIVE && find_set(uf, i) != root) {
            connected = 0;
            break;
        }
    }
    freeUnionFind(uf);
    return connected;
}

// --- MODIFIED: save_graph_to_file now uses the stored coordinates ---
void save_graph_to_file(Graph* g, const char* dims_str, long long trial) {
    char nodes_filename[256], edges_filename[256];
    // MODIFIED: Use %s for dims_str
    sprintf(nodes_filename, "output/network_N%s_trial%lld.nodes", dims_str, trial);
    sprintf(edges_filename, "output/network_N%s_trial%lld.edges", dims_str, trial);
    FILE* nodes_file = fopen(nodes_filename, "w");
    if (!nodes_file) { perror("Failed to open nodes file"); return; }
    /* Record the true periodic cell. Without it the Python loader has to
     * estimate the box from the coordinate extent as max-min+1, which is
     * exact only for SC: BCC/FCC basis sites sit at +0.5 and never reach
     * the cell edge, so the estimate overshoots by half a cell and sends
     * a large fraction of edges to the wrong periodic replica. Must stay
     * byte-compatible with topon.topology.loader.format_box_header. */
    {
        int bx = 0, by = 0, bz = 0;
        if (sscanf(dims_str, "%dx%dx%d", &bx, &by, &bz) == 3) {
            fprintf(nodes_file, "# BOX %g %g %g\n",
                    (double)bx, (double)by, (double)bz);
        }
    }
    fprintf(nodes_file, "# NodeID X Y Z Degree\n");
    for (int i = 0; i < g->V; ++i) {
        // Use the stored coordinates directly
        fprintf(nodes_file, "%d %f %f %f %d\n", i, g->coords[i].x, g->coords[i].y, g->coords[i].z, g->degrees[i]);
    }
    fclose(nodes_file);
    FILE* edges_file = fopen(edges_filename, "w");
    if (!edges_file) { perror("Failed to open edges file"); return; }
    fprintf(edges_file, "# Node1 Node2\n");
    for (int i = 0; i < g->V; ++i) {
        AdjListNode* pCrawl = g->array[i].head;
        while(pCrawl) {
            if (i < pCrawl->dest) fprintf(edges_file, "%d %d\n", i, pCrawl->dest);
            pCrawl = pCrawl->next;
        }
    }
    fclose(edges_file);
    printf("Successfully saved network from trial %lld to files.\n", trial);
}

void save_move_log_to_file(const MoveLog* move_log, long long count, const char* dims_str, long long trial) {
    char log_filename[256];
    // MODIFIED: Use %s for dims_str
    sprintf(log_filename, "output/network_N%s_trial%lld.log", dims_str, trial);
    FILE* log_file = fopen(log_filename, "w");
    if (!log_file) {
        perror("Failed to open move log file");
        return;
    }
    fprintf(log_file, "# Successful move log for Trial %lld\n", trial);
    for (long long i = 0; i < count; ++i) {
        const MoveLog* move = &move_log[i];
        switch (move->type) {
            case SET_D0:
                fprintf(log_file, "set node %d to d0\n", move->u);
                break;
            case SET_D1:
                fprintf(log_file, "set node %d to d1\n", move->u);
                break;
            case REMOVE_EDGE:
                fprintf(log_file, "remove edge %d-%d\n", move->u, move->v);
                break;
        }
    }
    fclose(log_file);
    printf("Successfully saved move log for trial %lld.\n", trial);
}


// --- Core Simulation Logic ---

void shuffle_array(int *array, size_t n) {
    if (n > 1) {
        for (size_t i = n - 1; i > 0; i--) {
            size_t j = rand() % (i + 1);
            int temp = array[i];
            array[i] = array[j];
            array[j] = temp;
        }
    }
}

// --- MODIFIED: Signature updated to parse 'e:N' ---
int parse_degree_distribution(char* str, int* target_counts, int max_degree_val, int* target_edge_count) {
    for(int i = 0; i <= max_degree_val; ++i) target_counts[i] = -2; // -2 means not specified
    char* token = strtok(str, ",");
    while (token != NULL) {
        int degree, count;
        // --- NEW: Check for e:N ---
        if (sscanf(token, "e:%d", &count) == 1) {
            *target_edge_count = count;
        }
        // --- END NEW ---
        else if (sscanf(token, "%d:%d", &degree, &count) == 2) {
            if (degree > max_degree_val) return 0; // Error
            target_counts[degree] = count;
        } else { return 0; } // Error
        token = strtok(NULL, ",");
    }
    return 1;
}

/**
 * @brief Prints the current degree distribution of the graph compared to the target.
 */
void print_distribution(const char* stage_name, Graph* g, const int* target_counts, long long trial_num, long long move_num, int max_func, int extensive_logging) {
    int max_current_degree = 0;
    for(int i=0; i < g->V; ++i) {
        if(g->degrees[i] > max_current_degree) max_current_degree = g->degrees[i];
    }
    
    int max_print_degree = max_current_degree > max_func ? max_current_degree : max_func;
    if (max_print_degree < 6) max_print_degree = 6;

    int buffer_size = max_print_degree + 1;
    int* current_counts = (int*)calloc(buffer_size, sizeof(int));
    long long current_total_degree_sum = 0; // --- NEW ---
    for(int i = 0; i < g->V; ++i) {
        if (g->degrees[i] < buffer_size) {
            current_counts[g->degrees[i]]++;
        }
        current_total_degree_sum += g->degrees[i]; // --- NEW ---
    }
    
    if (move_num > 0) {
         printf("[Trial %lld | Move %-8lld] ", trial_num, move_num);
    } else {
         printf("[Trial %lld | %-24s] ", trial_num, stage_name);
    }

    // --- NEW: Print total edge count ---
    printf("Edges: %-6lld | ", current_total_degree_sum / 2);
    printf("Dist:");
    for(int i = 0; i <= max_print_degree; ++i) {
        if (current_counts[i] > 0 || (i <= max_func && target_counts[i] != -2) ) {
            printf(" d%d:%d", i, current_counts[i]);
            if (i <= max_func && target_counts[i] != -2) {
                if (target_counts[i] == -1) printf("/*");
                else printf("/%d", target_counts[i]);
            }
        }
    }
    printf("\n");
    fflush(stdout);
    free(current_counts);
}
// --- MODIFIED: Added new 'lattice_type' and 'target_edge_count' arguments ---
Graph* run_single_trial(Graph* base_graph, int max_func, const int* target_counts, int target_edge_count, long long trial_num, int extensive_logging, const char* dims_str, const char* lattice_type) {
    Graph* g = copyGraph(base_graph);
    int total_nodes = g->V;
    NodeStatus* node_status = (NodeStatus*)malloc(total_nodes * sizeof(NodeStatus));
    int* node_indices = (int*)malloc(total_nodes * sizeof(int));
    Graph* return_graph = NULL;
    long long move_counter = 0;
    
    // --- NEW: Calculate target total degree sum from target edge count ---
    // -1 (from target_edge_count) * 2 = -2. This is our "not set" flag.
    long long target_degree_sum = (long long)target_edge_count * 2;

    // --- NEW: Check lattice type once at the beginning for efficiency ---
    int is_sc_lattice = (strcmp(lattice_type, "SC") == 0);

    MoveLog* move_log = NULL;
    long long move_log_count = 0;
    long long move_log_capacity = 0;
    if (extensive_logging >= 1) {
        move_log_capacity = g->V * 3; 
        move_log = (MoveLog*)malloc(move_log_capacity * sizeof(MoveLog));
    }

    for(int i=0; i<total_nodes; ++i) {
        node_status[i] = ACTIVE;
        node_indices[i] = i;
    }
    shuffle_array(node_indices, total_nodes);

    int N0_target = (target_counts[0] >= 0) ? target_counts[0] : 0;
    int N1_target = (target_counts[1] >= 0) ? target_counts[1] : 0;
    if (N0_target + N1_target > total_nodes) goto cleanup;
    
    int current_node_offset = 0;

    // --- Stage 1: Set Degree-0 Nodes (Strict) ---
    for(int i=0; i<N0_target; ++i) {
        int node_idx = node_indices[current_node_offset++];
        while(g->degrees[node_idx] > 0) {
            int num_neighbors = g->degrees[node_idx];
            int* neighbors = (int*)malloc(num_neighbors * sizeof(int));
            AdjListNode* pCrawl = g->array[node_idx].head;
            for(int k=0; k<num_neighbors; ++k){ neighbors[k] = pCrawl->dest; pCrawl = pCrawl->next; }
            
            int removed = 0;
            for(int k=0; k<num_neighbors; ++k) {
                int neighbor_idx = neighbors[k];

                // --- MODIFIED: This check is now conditional on being an SC lattice ---
                if (is_sc_lattice && g->degrees[neighbor_idx] <= 2) {
                    continue;
                }
                
                // --- MODIFIED: Pass sums to is_move_safe (not relevant for stage 1, pass -1)
                if (!is_move_safe(g, node_idx, neighbor_idx, target_counts, max_func, 1, target_degree_sum, -1)) {
                    continue;
                }
                
                if (extensive_logging >= 1) { /* log move */ }
                removeEdge(g, node_idx, neighbor_idx);
                
                if (extensive_logging == 1) { /* print per-move */ }
                removed = 1;
                break;
            }
            free(neighbors);
            if(!removed) goto cleanup; 
        }
        node_status[node_idx] = IS_DEGREE_0;
        if (extensive_logging >= 1) { /* log move */ }
    }
    if (extensive_logging >= 1) print_distribution("Stage 1: Set-d0 (Strict)", g, target_counts, trial_num, 0, max_func, extensive_logging);

    // --- Stage 2: Set Degree-1 Nodes (Strict) ---
    for(int i=0; i<N1_target; ++i) {
       int node_idx = node_indices[current_node_offset++];
       if(node_status[node_idx] != ACTIVE) { i--; continue; } 
       
       while(g->degrees[node_idx] > 1) {
            int num_neighbors = g->degrees[node_idx];
            int* neighbors = (int*)malloc(num_neighbors * sizeof(int));
            AdjListNode* pCrawl = g->array[node_idx].head;
            for(int k=0; k<num_neighbors; ++k){ neighbors[k] = pCrawl->dest; pCrawl = pCrawl->next; }
            shuffle_array(neighbors, num_neighbors);
            
            int removed = 0;
            for(int k=0; k<num_neighbors; ++k){
                int neighbor_idx = neighbors[k];
                
                // --- MODIFIED: This check is now conditional on being an SC lattice ---
                if (is_sc_lattice && g->degrees[neighbor_idx] <= 2) {
                    continue;
                }

                // --- MODIFIED: Pass sums to is_move_safe (not relevant for stage 2, pass -1)
                if (!is_move_safe(g, node_idx, neighbor_idx, target_counts, max_func, 2, target_degree_sum, -1)) {
                    continue;
                }
                
                if (extensive_logging >= 1) { /* log move */ }
                removeEdge(g, node_idx, neighbor_idx);

                if(is_subgraph_connected(g, node_status)){
                    if (extensive_logging == 1) { /* print per-move */ }
                    removed = 1;
                    break;
                }
                else { 
                    addEdge(g, node_idx, neighbor_idx); // Backtrack
                    if(extensive_logging >= 1) move_log_count--;
                } 
            }
            free(neighbors);
            if(!removed) goto cleanup;
       }
       node_status[node_idx] = IS_DEGREE_1;
       if (extensive_logging >= 1) { /* log move */ }
    }
    if (extensive_logging >= 1) print_distribution("Stage 2: Set-d1 (Strict)", g, target_counts, trial_num, 0, max_func, extensive_logging);

    // --- Stage 3: Enforce Max Functionality (Strict) ---
    for(int i=0; i<total_nodes; ++i) {
        int node_idx = node_indices[i]; 
        if(node_status[node_idx] != ACTIVE) continue; 
        
        while(g->degrees[node_idx] > max_func) {
            int num_neighbors = g->degrees[node_idx];
            int* neighbors = (int*)malloc(num_neighbors * sizeof(int));
            AdjListNode* pCrawl = g->array[node_idx].head;
            for(int k=0; k<num_neighbors; ++k){ neighbors[k] = pCrawl->dest; pCrawl = pCrawl->next; }
            shuffle_array(neighbors, num_neighbors);

            int removed = 0;
            for(int k=0; k<num_neighbors; ++k) {
                int neighbor_idx = neighbors[k];

                // --- MODIFIED: This check is now conditional on being an SC lattice ---
                if (is_sc_lattice && g->degrees[neighbor_idx] <= 2) {
                    continue;
                }

                // --- MODIFIED: Pass sums to is_move_safe (not relevant for stage 3, pass -1)
                if (!is_move_safe(g, node_idx, neighbor_idx, target_counts, max_func, 3, target_degree_sum, -1)) {
                    continue;
                }

                if (extensive_logging >= 1) {
                    if (move_log_count >= move_log_capacity) {
                        move_log_capacity *= 2;
                        move_log = (MoveLog*)realloc(move_log, move_log_capacity * sizeof(MoveLog));
                    }
                    move_log[move_log_count++] = (MoveLog){REMOVE_EDGE, node_idx, neighbor_idx};
                }
                removeEdge(g, node_idx, neighbor_idx);

                if (is_subgraph_connected(g, node_status)) {
                    if (extensive_logging == 1) {
                        move_counter++;
                        print_distribution("Stage 3: Enforce-Max", g, target_counts, trial_num, move_counter, max_func, extensive_logging);
                    }
                    removed = 1;
                    break;
                } else {
                    addEdge(g, node_idx, neighbor_idx);
                    if(extensive_logging >= 1) move_log_count--;
                }
            }
            free(neighbors);
            if (!removed) goto cleanup;
        }
    }
    if (extensive_logging >= 1) print_distribution("Stage 3: Enforce-Max (Strict)", g, target_counts, trial_num, 0, max_func, extensive_logging);

    // --- Stage 4: Systematic Search Loop ---
    while(1) {
        int* current_counts = (int*)calloc(max_func + 3, sizeof(int));
        int is_done = 1;
        long long current_total_degree_sum = 0; // --- NEW ---
        
        int has_high_degree_nodes = 0;
        for(int i=0; i<total_nodes; ++i) {
             current_total_degree_sum += g->degrees[i]; // --- NEW ---
             if (g->degrees[i] <= max_func + 2) current_counts[g->degrees[i]]++;
             if (node_status[i] == ACTIVE && g->degrees[i] > max_func) {
                 has_high_degree_nodes = 1;
             }
        }

        // ---
        // --- THIS IS THE CORRECTED LOGIC ---
        // ---
        if (has_high_degree_nodes) {
            is_done = 0; // Not done, still have nodes > max_func (e.g., in mf=4 run)
        } else {
            // 1. Check all *explicit* d:N targets first
            for(int i=0; i <= max_func; ++i) {
                if (target_counts[i] >= 0 && target_counts[i] != current_counts[i]) {
                    is_done = 0; // Failed an explicit target
                    break;
                }
            }

            if (is_done) { 
                // 2. If explicit targets are met, check which mode we are in
                if (target_edge_count != -1) {
                    // --- e:N Mode (Bond Percolation) ---
                    // This is for your bond percolation study.
                    if (current_total_degree_sum != target_degree_sum) {
                        is_done = 0; // Edge count is wrong
                    }
                    // We also must be connected (implicitly checked by Stage 4 moves)
                    if (!is_subgraph_connected(g, node_status)) {
                         is_done = 0; // Not connected
                    }
                } else {
                    // --- Legacy Mode (Site Percolation) ---
                    // This is for your site percolation study.
                    // We have already confirmed explicit targets (d0) are met.
                    // We *allow* d1-d6 to be "anything".
                    // The *only* thing to check now is connectivity.
                    
                    if (is_subgraph_connected(g, node_status)) {
                        // SUCCESS! Explicit targets (d0) met AND it's connected.
                        is_done = 1;
                    } else {
                        // FAILED! Explicit targets (d0) met but NOT connected.
                        // This is a "valid" data point, but not a successful trial.
                        free(current_counts);
                        goto cleanup; // Fail the trial
                    }
                }
            }
        }
        // --- END CORRECTED LOGIC ---
        
        if (is_done) {
            free(current_counts);
            print_distribution("Final Distribution", g, target_counts, trial_num, 0, max_func, extensive_logging);
            return_graph = g;
            if (extensive_logging >= 1) {
                save_move_log_to_file(move_log, move_log_count, dims_str, trial_num);
            }
            goto cleanup_success;
        }
        free(current_counts);

        // ---
        // This part of the loop will now only be reached if:
        // 1. We are in 'e:N' mode and haven't hit the target edge count.
        // 2. We are in 'mf=4' mode and haven't finished Stage 3's job.
        // It will *not* be reached by a 'mf=6' site percolation run
        // because that run will either 'goto cleanup' or 'goto cleanup_success'.
        // ---
        if (current_total_degree_sum == 0) goto cleanup;
        
        long long num_edges = current_total_degree_sum / 2;
        Edge* edge_list = (Edge*)malloc(num_edges * sizeof(Edge));
        long long current_edge_idx = 0;
        // --- Later we can add spatial skew?
        for(int j=0; j<total_nodes; ++j) { 
            int i = node_indices[j]; // Use the shuffled index
            if (node_status[i] == ACTIVE) {
                AdjListNode* pCrawl = g->array[i].head;
                while(pCrawl) {
                    // We can still use i < pCrawl->dest to avoid double counting
                    if (i < pCrawl->dest && node_status[pCrawl->dest] == ACTIVE) { 
                        edge_list[current_edge_idx].u = i;
                        edge_list[current_edge_idx].v = pCrawl->dest;
                        current_edge_idx++;
                    }
                    pCrawl = pCrawl->next;
                }
            }
        }
        shuffle_array((int*)edge_list, current_edge_idx);

        int move_made = 0;
        for(long long i=0; i<current_edge_idx; ++i) {
            int u = edge_list[i].u;
            int v = edge_list[i].v;

            int u_deg = g->degrees[u];
            int v_deg = g->degrees[v];

            if (u_deg <= 1 || v_deg <= 1) continue;
            
            // --- MODIFIED: Legacy check for d2, still useful ---
            if (u_deg == 2 || v_deg == 2) {
                // This logic is for legacy mode, but is also safe for 'e:N' mode.
                // It prevents over-sculpting d1 if d1 is not an explicit target.
                // If d1 target is 0, is_move_safe will catch it anyway.
                if (target_counts[1] != -1) { 
                    int current_d1_count = 0;
                    for (int j = 0; j < total_nodes; j++) {
                        if (node_status[j] == ACTIVE && g->degrees[j] == 1) current_d1_count++;
                    }
                    if (current_d1_count >= target_counts[1]) {
                        continue; 
                    }
                }
            }
            
            // --- MODIFIED: Pass current degree sum to is_move_safe ---
            if (!is_move_safe(g, u, v, target_counts, max_func, 4, target_degree_sum, current_total_degree_sum)) {
                continue;
            }

            if (extensive_logging >= 1) { /* log move */ }
            removeEdge(g, u, v);
            
            if (is_subgraph_connected(g, node_status)) {
                move_made = 1;
                break;
            } else {
                addEdge(g, u, v);
                if(extensive_logging >= 1) move_log_count--;
            }
        }
        free(edge_list);

        if (move_made) {
            move_counter++;
            if (extensive_logging == 1) {
                print_distribution("Stage 4: Systematic Search", g, target_counts, trial_num, move_counter, max_func, extensive_logging);
            }
        } else {
            // --- NEW: If no moves are safe, but we are in e:N mode, check if we are simply stuck
            if (target_edge_count != -1) {
                // We are stuck. We already know 'is_done' is false from the start of the loop.
                // This means we are either stuck *above* the target edge count (failure)
                // or *at* the target edge count but with wrong d:N explicit targets (failure).
                goto cleanup;
            }
            
            // Legacy mode: no moves means failure
            goto cleanup;
        }
    }

cleanup:
    freeGraph(g);
    g = NULL;
cleanup_success:
    if (extensive_logging >= 1) { free(move_log); }
    free(node_status);
    free(node_indices);
    return g;
}

// --- START: LATTICE CREATION LOGIC ---

// Helper for SC lattice (original logic)
Graph* create_sc_lattice(int Nx, int Ny, int Nz, const int* p_dims) {
    int total_nodes = Nx * Ny * Nz; 
    Graph* g = createGraph(total_nodes);
    
    // C++ Lambda REMOVED: auto get_index = [&](int x, int y, int z) { ... };

    for (int z = 0; z < Nz; z++) { 
        for (int y = 0; y < Ny; y++) { 
            for (int x = 0; x < Nx; x++) { 
                int u = x + y * Nx + z * Nx * Ny; 
                g->coords[u] = (Coord){(double)x, (double)y, (double)z};

                if (p_dims[0] || x < Nx - 1) addEdge(g, u, ((x + 1) % Nx) + y * Nx + z * Nx * Ny);
                if (p_dims[1] || y < Ny - 1) addEdge(g, u, x + ((y + 1) % Ny) * Nx + z * Nx * Ny);
                if (p_dims[2] || z < Nz - 1) addEdge(g, u, x + y * Nx + ((z + 1) % Nz) * Nx * Ny);
            }
        }
    }
    return g;
}

// Helper for BCC lattice
Graph* create_bcc_lattice(int Nx, int Ny, int Nz, const int* p_dims) {
    int total_nodes = 2 * Nx * Ny * Nz; // MODIFIED
    Graph* g = createGraph(total_nodes);
    int node_idx = 0;

    int high_res_Nx = 2 * Nx; // MODIFIED
    int high_res_Ny = 2 * Ny; // MODIFIED
    int high_res_Nz = 2 * Nz; // MODIFIED
    // MODIFIED: Use new high_res dims
    long long map_size = (long long)high_res_Nx * high_res_Ny * high_res_Nz;
    int* coord_to_id_map = (int*)malloc(map_size * sizeof(int));
    for(long long i=0; i<map_size; ++i) coord_to_id_map[i] = -1;

    // C++ Lambda REMOVED: auto get_map_idx = [&](int x, int y, int z) { ... };

    // 1. Place nodes
    for (int k = 0; k < Nz; k++) { // MODIFIED
        for (int j = 0; j < Ny; j++) { // MODIFIED
            for (int i = 0; i < Nx; i++) { // MODIFIED
                // Corner node
                int cx = 2*i, cy = 2*j, cz = 2*k;
                g->coords[node_idx] = (Coord){(double)i, (double)j, (double)k};
                // C equivalent: use direct calculation (MODIFIED for new indexing)
                coord_to_id_map[(long long)cx + (long long)cy * high_res_Nx + (long long)cz * high_res_Nx * high_res_Ny] = node_idx++;
                
                // Body-centered node
                int bx = 2*i+1, by = 2*j+1, bz = 2*k+1;
                g->coords[node_idx] = (Coord){(double)i+0.5, (double)j+0.5, (double)k+0.5};
                // C equivalent: use direct calculation (MODIFIED for new indexing)
                coord_to_id_map[(long long)bx + (long long)by * high_res_Nx + (long long)bz * high_res_Nx * high_res_Ny] = node_idx++;
            }
        }
    }

    // 2. Connect nodes (8 nearest neighbors)
    for (long long map_idx = 0; map_idx < map_size; ++map_idx) {
        int id = coord_to_id_map[map_idx];
        if (id == -1) continue;

        // MODIFIED: Use new high_res dims for coordinate extraction
        int z = map_idx / ((long long)high_res_Nx * high_res_Ny);
        int y = (map_idx / high_res_Nx) % high_res_Ny;
        int x = map_idx % high_res_Nx;

        for (int dz = -1; dz <= 1; dz += 2) {
            for (int dy = -1; dy <= 1; dy += 2) {
                for (int dx = -1; dx <= 1; dx += 2) {
                    int nx = x + dx;
                    int ny = y + dy;
                    int nz = z + dz;
                    // Handle periodicity (MODIFIED for new high_res dims)
                    if (p_dims[0]) nx = (nx + high_res_Nx) % high_res_Nx;
                    if (p_dims[1]) ny = (ny + high_res_Ny) % high_res_Ny;
                    if (p_dims[2]) nz = (nz + high_res_Nz) % high_res_Nz;

                    // MODIFIED for new high_res dims
                    if (nx >= 0 && nx < high_res_Nx && ny >= 0 && ny < high_res_Ny && nz >= 0 && nz < high_res_Nz) {
                        // C equivalent: use direct calculation (MODIFIED for new indexing)
                        int neighbor_id = coord_to_id_map[(long long)nx + (long long)ny * high_res_Nx + (long long)nz * high_res_Nx * high_res_Ny];
                        if (neighbor_id != -1 && id < neighbor_id) {
                            addEdge(g, id, neighbor_id);
                        }
                    }
                }
            }
        }
    }
    free(coord_to_id_map);
    return g;
}

// Helper for FCC lattice
Graph* create_fcc_lattice(int Nx, int Ny, int Nz, const int* p_dims) {
    int total_nodes = 4 * Nx * Ny * Nz; // MODIFIED
    Graph* g = createGraph(total_nodes);
    int node_idx = 0;

    int high_res_Nx = 2 * Nx; // MODIFIED
    int high_res_Ny = 2 * Ny; // MODIFIED
    int high_res_Nz = 2 * Nz; // MODIFIED
    // MODIFIED: Use new high_res dims
    long long map_size = (long long)high_res_Nx * high_res_Ny * high_res_Nz;
    int* coord_to_id_map = (int*)malloc(map_size * sizeof(int));
    for(long long i=0; i<map_size; ++i) coord_to_id_map[i] = -1;


    // 1. Place nodes
    for (int k = 0; k < Nz; k++) { // MODIFIED
        for (int j = 0; j < Ny; j++) { // MODIFIED
            for (int i = 0; i < Nx; i++) { // MODIFIED
                // Corner node
                g->coords[node_idx] = (Coord){(double)i, (double)j, (double)k};
                // MODIFIED: New indexing
                coord_to_id_map[(long long)(2*i) + (long long)(2*j) * high_res_Nx + (long long)(2*k) * high_res_Nx * high_res_Ny] = node_idx++;
                // Face nodes
                g->coords[node_idx] = (Coord){(double)i+0.5, (double)j+0.5, (double)k};
                // MODIFIED: New indexing
                coord_to_id_map[(long long)(2*i+1) + (long long)(2*j+1) * high_res_Nx + (long long)(2*k) * high_res_Nx * high_res_Ny] = node_idx++;
                g->coords[node_idx] = (Coord){(double)i+0.5, (double)j, (double)k+0.5};
                // MODIFIED: New indexing
                coord_to_id_map[(long long)(2*i+1) + (long long)(2*j) * high_res_Nx + (long long)(2*k+1) * high_res_Nx * high_res_Ny] = node_idx++;
                g->coords[node_idx] = (Coord){(double)i, (double)j+0.5, (double)k+0.5};
                // MODIFIED: New indexing
                coord_to_id_map[(long long)(2*i) + (long long)(2*j+1) * high_res_Nx + (long long)(2*k+1) * high_res_Nx * high_res_Ny] = node_idx++;
            }
        }
    }

// 2. Connect nodes (12 nearest neighbors)
    for (long long map_idx = 0; map_idx < map_size; ++map_idx) {
        int id = coord_to_id_map[map_idx];
        if (id == -1) continue;

        // MODIFIED: Use new high_res dims for coordinate extraction
        int z = map_idx / ((long long)high_res_Nx * high_res_Ny);
        int y = (map_idx / high_res_Nx) % high_res_Ny;
        int x = map_idx % high_res_Nx;

        // --- MODIFIED: Define all 12 neighbor directions ---
        // (No change to this array)
        int neighbor_offsets[12][3] = {
            {1,1,0}, {1,-1,0}, {-1,1,0}, {-1,-1,0},  // XY plane
            {1,0,1}, {1,0,-1}, {-1,0,1}, {-1,0,-1},  // XZ plane
            {0,1,1}, {0,1,-1}, {0,-1,1}, {0,-1,-1}   // YZ plane
        };

        // --- MODIFIED: Loop over all 12 offsets ---
        for(int i=0; i<12; ++i) {
            int nx = x + neighbor_offsets[i][0];
            int ny = y + neighbor_offsets[i][1];
            int nz = z + neighbor_offsets[i][2];
            
            // Handle periodicity (MODIFIED for new high_res dims)
            if (p_dims[0]) nx = (nx + high_res_Nx) % high_res_Nx;
            if (p_dims[1]) ny = (ny + high_res_Ny) % high_res_Ny;
            if (p_dims[2]) nz = (nz + high_res_Nz) % high_res_Nz;

            // MODIFIED for new high_res dims
            if (nx >= 0 && nx < high_res_Nx && ny >= 0 && ny < high_res_Ny && nz >= 0 && nz < high_res_Nz) {
                // MODIFIED: New indexing
                int neighbor_id = coord_to_id_map[(long long)nx + (long long)ny * high_res_Nx + (long long)nz * high_res_Nx * high_res_Ny];
                
                // The (id < neighbor_id) check correctly prevents double-counting
                if (neighbor_id != -1 && id < neighbor_id) { 
                    addEdge(g, id, neighbor_id);
                }
            }
        }
    }
    free(coord_to_id_map);
    return g;
}

/* Diamond lattice: 8 sites per conventional cubic cell, every site
 * exactly 4-coordinated by construction. Two interpenetrating FCC
 * sublattices offset by (1/4, 1/4, 1/4) along the body diagonal.
 *
 * This is the cleanest backbone for a max_func=4 network: the raw
 * lattice already satisfies that ceiling, so sculpting has nothing to
 * prune unless defects are requested.
 *
 * Mirrors create_diamond_lattice in generator_python_diamond.py,
 * including the site order (cell-major, then the eight basis sites in
 * the order listed below), so the two number their nodes identically.
 *
 * Sites live on a x4 integer grid. A-sublattice sites satisfy
 * (hx+hy+hz) % 4 == 0 and reach their four neighbours through the
 * sign-product +1 offsets; B-sublattice sites satisfy sum % 4 == 3 and
 * use the sign-product -1 offsets. No other residue holds a site.
 */
Graph* create_diamond_lattice(int Nx, int Ny, int Nz, const int* p_dims) {
    static const int basis_hr[8][3] = {
        {0,0,0}, {2,2,0}, {2,0,2}, {0,2,2},     /* A sublattice */
        {1,1,1}, {3,3,1}, {3,1,3}, {1,3,3},     /* B sublattice */
    };
    static const int off_a[4][3] = {            /* sign product +1: A -> B */
        {+1,+1,+1}, {+1,-1,-1}, {-1,+1,-1}, {-1,-1,+1},
    };
    static const int off_b[4][3] = {            /* sign product -1: B -> A */
        {-1,-1,-1}, {-1,+1,+1}, {+1,-1,+1}, {+1,+1,-1},
    };

    int total_nodes = 8 * Nx * Ny * Nz;
    Graph* g = createGraph(total_nodes);
    int node_idx = 0;

    int hr_Nx = 4 * Nx, hr_Ny = 4 * Ny, hr_Nz = 4 * Nz;
    long long map_size = (long long)hr_Nx * hr_Ny * hr_Nz;
    int* coord_to_id_map = (int*)malloc(map_size * sizeof(int));
    if (!coord_to_id_map) {
        fprintf(stderr, "Error: out of memory building diamond lattice.\n");
        freeGraph(g);
        return NULL;
    }
    for (long long i = 0; i < map_size; ++i) coord_to_id_map[i] = -1;

    /* 1. Place sites. */
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                for (int b = 0; b < 8; ++b) {
                    int hx = 4*i + basis_hr[b][0];
                    int hy = 4*j + basis_hr[b][1];
                    int hz = 4*k + basis_hr[b][2];
                    g->coords[node_idx] = (Coord){
                        (double)i + basis_hr[b][0] / 4.0,
                        (double)j + basis_hr[b][1] / 4.0,
                        (double)k + basis_hr[b][2] / 4.0,
                    };
                    coord_to_id_map[(long long)hx
                                    + (long long)hy * hr_Nx
                                    + (long long)hz * hr_Nx * hr_Ny] = node_idx++;
                }
            }
        }
    }

    /* 2. Connect each site to its four tetrahedral neighbours. */
    for (long long map_idx = 0; map_idx < map_size; ++map_idx) {
        int id = coord_to_id_map[map_idx];
        if (id == -1) continue;

        int z = map_idx / ((long long)hr_Nx * hr_Ny);
        int y = (map_idx / hr_Nx) % hr_Ny;
        int x = map_idx % hr_Nx;

        int s = (x + y + z) % 4;
        const int (*offsets)[3] = (s == 0) ? off_a : off_b;

        for (int o = 0; o < 4; ++o) {
            int nx = x + offsets[o][0];
            int ny = y + offsets[o][1];
            int nz = z + offsets[o][2];
            if (p_dims[0]) nx = (nx + hr_Nx) % hr_Nx;
            if (p_dims[1]) ny = (ny + hr_Ny) % hr_Ny;
            if (p_dims[2]) nz = (nz + hr_Nz) % hr_Nz;

            if (nx >= 0 && nx < hr_Nx && ny >= 0 && ny < hr_Ny && nz >= 0 && nz < hr_Nz) {
                int neighbor_id = coord_to_id_map[(long long)nx
                                                  + (long long)ny * hr_Nx
                                                  + (long long)nz * hr_Nx * hr_Ny];
                if (neighbor_id != -1 && id < neighbor_id) {
                    addEdge(g, id, neighbor_id);
                }
            }
        }
    }

    free(coord_to_id_map);
    return g;
}

/* Mixed SC/BCC/FCC lattice.
 *
 * Mirrors PythonTopologyGenerator._create_mixed_lattice. All three
 * lattices share the cubic cell corner and each adds sites on top of it:
 * BCC one body centre, FCC three face centres. So the corner goes into
 * every cell, the body centre with probability f_bcc and each face
 * centre with probability f_fcc. The SC fraction is the remainder and
 * contributes no site of its own, which is what makes the three
 * fractions a partition summing to 1.
 *
 * Edges join every pair within `cutoff` under the minimum image, because
 * a mixed point set has no single neighbour shell to enumerate the way
 * the pure builders do. Non-periodic axes (p_dims[k] == 0) are not
 * wrapped. The search is O(N^2) in the site count, which is negligible
 * at the cell counts topon uses but grows quickly past ~20x20x20.
 *
 * Site order matches the Python builder (z outer, then y, then x, and
 * within a cell corner, body, then the XY/XZ/YZ faces) so that fractions
 * (1, 0, 0) reproduce create_sc_lattice exactly. The draws come from
 * rand(), so C and Python agree on distributions, not on a given draw.
 */
Graph* create_mixed_lattice(int Nx, int Ny, int Nz, const int* p_dims,
                            double f_bcc, double f_fcc, double cutoff) {
    static const double face_off[3][3] = {
        {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}
    };

    /* Worst case is every optional site drawn: corner + body + 3 faces = 5
     * per cell. Sizing at 4 (the FCC site count) is NOT enough: with
     * f_bcc=0.1, f_fcc=0.9 a 4x4x4 overruns 4*N in ~0.2% of draws. */
    int max_sites = 5 * Nx * Ny * Nz;
    Coord* sites = (Coord*)malloc((size_t)max_sites * sizeof(Coord));
    if (!sites) { fprintf(stderr, "Error: out of memory building mixed lattice.\n"); return NULL; }

    int n = 0;
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                sites[n++] = (Coord){(double)i, (double)j, (double)k};
                if (f_bcc > 0.0 && ((double)rand() / ((double)RAND_MAX + 1.0)) < f_bcc) {
                    sites[n++] = (Coord){i + 0.5, j + 0.5, k + 0.5};
                }
                if (f_fcc > 0.0) {
                    for (int f = 0; f < 3; ++f) {
                        if (((double)rand() / ((double)RAND_MAX + 1.0)) < f_fcc) {
                            sites[n++] = (Coord){i + face_off[f][0],
                                                 j + face_off[f][1],
                                                 k + face_off[f][2]};
                        }
                    }
                }
            }
        }
    }

    Graph* g = createGraph(n);
    for (int i = 0; i < n; ++i) g->coords[i] = sites[i];

    double box[3] = {(double)Nx, (double)Ny, (double)Nz};
    double cutoff_sq = cutoff * cutoff;
    for (int a = 0; a < n; ++a) {
        for (int b = a + 1; b < n; ++b) {
            double d2 = 0.0;
            double dv[3] = {sites[a].x - sites[b].x,
                            sites[a].y - sites[b].y,
                            sites[a].z - sites[b].z};
            for (int ax = 0; ax < 3; ++ax) {
                double d = dv[ax];
                if (p_dims[ax]) {
                    /* Minimum image; rint() rounds half to even, matching
                     * numpy.round in the Python builder. */
                    d -= box[ax] * rint(d / box[ax]);
                }
                d2 += d * d;
            }
            if (d2 <= cutoff_sq + 1e-12 && d2 > 1e-12) addEdge(g, a, b);
        }
    }

    free(sites);
    return g;
}

// --- END: LATTICE CREATION LOGIC ---


// --- Main ---

int main(int argc, char *argv[]) {
    // --- MODIFIED: Check for 9 arguments ---
    if (argc != 9) {
        // MODIFIED: Updated usage string
        fprintf(stderr, "Usage: %s <dims_str> <periodicity_str> <max_func> <max_trials> <max_saves> \"<degree_dist_string>\" <extensive_logging> <lattice_type>\n", argv[0]);
        fprintf(stderr, "Example (Legacy): %s 8x8x8 111 4 1000 1 \"0:0,1:0,2:100,3:312,4:100\" 1 FCC\n", argv[0]);
        fprintf(stderr, "Example (New 'e'): %s 8x6x8 110 6 1000 1 \"0:0,1:0,2:100,e:450\" 1 SC\n", argv[0]);
        fprintf(stderr, "  <dims_str>: Dimensions in NxN_yN_z format (e.g., '8x6x8').\n");
        fprintf(stderr, "  <degree_dist_string>: \"d0:N0,d1:N1,e:TotalEdges\"\n");
        fprintf(stderr, "  <lattice_type>: SC, BCC, FCC, Diamond, or MIX:<sc>,<bcc>,<fcc>[,<cutoff>]\n");
        fprintf(stderr, "Example (Mixed):  %s 6x6x6 111 4 1000 1 \"0:0,1:0\" 0 MIX:0.2,0.4,0.4\n", argv[0]);
        return 1;
    }

    // MODIFIED: Parse dims_str instead of N
    char* dims_str = argv[1];
    int Nx, Ny, Nz;
    if (sscanf(dims_str, "%dx%dx%d", &Nx, &Ny, &Nz) != 3) {
        fprintf(stderr, "Error: Invalid dimensions string '%s'. Expected format: NxN_yN_z (e.g., '8x8x8' or '8x6x8').\n", dims_str);
        return 1;
    }

    char* periodicity_str = argv[2];
    int max_func = atoi(argv[3]);
    int max_trials = atoi(argv[4]);
    int max_saves = atoi(argv[5]);
    char* degree_dist_string_arg = argv[6];
    int extensive_logging = atoi(argv[7]);
    // --- NEW: Parse lattice type ---
    char* lattice_type = argv[8];

    int p_dims[3];
    p_dims[0] = (periodicity_str[0] == '1');
    p_dims[1] = (periodicity_str[1] == '1');
    p_dims[2] = (periodicity_str[2] == '1');
    
    char degree_dist_string[1024];
    strncpy(degree_dist_string, degree_dist_string_arg, sizeof(degree_dist_string) - 1);
    degree_dist_string[sizeof(degree_dist_string) - 1] = '\0';

    printf("SIMULATION: Starting serial execution with Strict Sculpting Algorithm.\n");
    // MODIFIED: Updated info print
    printf("INFO: Dims=%s, periodicity=[%d,%d,%d], max_func=%d, trials=%d, max_saves=%d, logging=%d, lattice=%s\n",
           dims_str, p_dims[0], p_dims[1], p_dims[2], max_func, max_trials, max_saves, extensive_logging, lattice_type);
    mkdir("output", 0755);

    /* Seed before building the lattice: MIX draws its sites from rand(),
     * unlike the pure builders. Moving this up from just before the trial
     * loop is a no-op for SC/BCC/FCC, which consume no randomness while
     * being built, so the sculpting stream they see is unchanged.
     *
     * The pid is mixed in because time(NULL) only advances once a second:
     * a script looping this executable to collect N networks used to get
     * byte-identical output from every run that started within the same
     * second. Verified before the fix -- three back-to-back runs produced
     * the same file. TOPON_SEED overrides for reproducible runs, which
     * the Python generator gets from seeding `random` directly. */
    unsigned int seed;
    const char* seed_env = getenv("TOPON_SEED");
    if (seed_env && *seed_env) {
        seed = (unsigned int)strtoul(seed_env, NULL, 10);
        printf("INFO: Using TOPON_SEED=%u (reproducible run).\n", seed);
    } else {
        seed = (unsigned int)time(NULL) ^ ((unsigned int)topon_getpid() << 16);
    }
    srand(seed);

    // --- MODIFIED: Create base_graph based on lattice_type input using Nx, Ny, Nz ---
    Graph* base_graph = NULL;
    if (strcmp(lattice_type, "SC") == 0) {
        base_graph = create_sc_lattice(Nx, Ny, Nz, p_dims);
    } else if (strcmp(lattice_type, "BCC") == 0) {
        base_graph = create_bcc_lattice(Nx, Ny, Nz, p_dims);
    } else if (strcmp(lattice_type, "FCC") == 0) {
        base_graph = create_fcc_lattice(Nx, Ny, Nz, p_dims);
    } else if (strcmp(lattice_type, "Diamond") == 0 || strcmp(lattice_type, "DIAMOND") == 0) {
        /* Both spellings: "Diamond" is what generator_python_diamond.py
         * uses, "DIAMOND" matches the shouting style of the other three. */
        base_graph = create_diamond_lattice(Nx, Ny, Nz, p_dims);
        if (!base_graph) return 1;
    } else if (strncmp(lattice_type, "MIX", 3) == 0 &&
               (lattice_type[3] == '\0' || lattice_type[3] == ':')) {
        /* The trailing check matters: a bare strncmp would also swallow
         * "MIXED", "MIXTURE" and any other typo starting with MIX, and
         * silently build a pure-SC lattice instead of rejecting it. */
        /* "MIX:<sc>,<bcc>,<fcc>[,<cutoff>]" -- the fractions ride inside
         * argv[8] so the positional CLI stays the shape every existing
         * caller and SLURM script already writes. */
        double f_sc = 1.0, f_bcc = 0.0, f_fcc = 0.0, cutoff = 1.0;
        const char* spec = lattice_type + 3;
        if (*spec == ':') {
            int got = sscanf(spec + 1, "%lf,%lf,%lf,%lf",
                             &f_sc, &f_bcc, &f_fcc, &cutoff);
            if (got < 3) {
                fprintf(stderr, "Error: MIX needs three fractions, e.g. MIX:0.2,0.4,0.4 "
                                "(optionally MIX:0.2,0.4,0.4,1.0 to set the cutoff). Got '%s'.\n",
                        lattice_type);
                return 1;
            }
        }
        double total = f_sc + f_bcc + f_fcc;
        if (f_sc < 0.0 || f_bcc < 0.0 || f_fcc < 0.0) {
            fprintf(stderr, "Error: MIX fractions must be non-negative (got %g,%g,%g).\n",
                    f_sc, f_bcc, f_fcc);
            return 1;
        }
        if (fabs(total - 1.0) > 1e-6) {
            fprintf(stderr, "Error: MIX fractions must sum to 1, got %g from %g,%g,%g. "
                            "They partition the crosslinker population, so a sum below 1 "
                            "would thin the lattice and above 1 would over-fill it.\n",
                    total, f_sc, f_bcc, f_fcc);
            return 1;
        }
        if (cutoff <= 0.0) {
            fprintf(stderr, "Error: MIX cutoff must be positive (got %g).\n", cutoff);
            return 1;
        }
        printf("INFO: Mixed lattice fractions SC=%g BCC=%g FCC=%g, cutoff=%g\n",
               f_sc, f_bcc, f_fcc, cutoff);
        base_graph = create_mixed_lattice(Nx, Ny, Nz, p_dims, f_bcc, f_fcc, cutoff);
        if (!base_graph) return 1;
    } else {
        fprintf(stderr, "Error: Invalid lattice type '%s'. Must be SC, BCC, FCC, Diamond, or MIX:<sc>,<bcc>,<fcc>.\n", lattice_type);
        return 1;
    }

    /* Size target_counts from the graph that was actually built, not from
     * a per-lattice constant. run_single_trial indexes this array by node
     * degree (target_counts[v_new_degree]), so a ceiling below the real
     * maximum is an out-of-bounds read. The pure lattices top out at a
     * known 6/8/12, but a MIX at the default cutoff already reaches 20,
     * and mix_cutoff is user-settable, so no constant is safe. The floor
     * of 12 preserves the array size the pure lattices always had. */
    int base_max_degree = 0;
    for (int i = 0; i < base_graph->V; ++i) {
        if (base_graph->degrees[i] > base_max_degree) base_max_degree = base_graph->degrees[i];
    }
    int max_possible_degree = max_func;
    if (base_max_degree > max_possible_degree) max_possible_degree = base_max_degree;
    if (max_possible_degree < 12) max_possible_degree = 12;

    int* target_counts = (int*)malloc((max_possible_degree + 1) * sizeof(int));
    int target_edge_count = -1; // --- NEW: For e:N ---

    char* str_copy = strdup(degree_dist_string);
    // --- MODIFIED: Pass target_edge_count pointer
    if (!parse_degree_distribution(str_copy, target_counts, max_possible_degree, &target_edge_count)) {
         fprintf(stderr, "Error parsing degree distribution string.\n");
         free(str_copy);
         free(target_counts);
         freeGraph(base_graph);
         return 1;
    }
    free(str_copy);

    long long node_sum = 0, explicit_degree_sum = 0; // --- MODIFIED: Renamed degree_sum
    int has_unconstrained = 0;
    for(int i=0; i <= max_func; ++i) {
        if(target_counts[i] >= 0) node_sum += target_counts[i];
        if(target_counts[i] == -1) has_unconstrained = 1;
        if(target_counts[i] > 0) explicit_degree_sum += (long long)i * target_counts[i];
    }

    /* Reject explicit targets for degrees above max_func.
     *
     * A node's final degree can never exceed max_func: stage 3 prunes to
     * it, and stage 4 refuses to finish while any ACTIVE node sits above
     * it. But the completion check only scans i <= max_func, so a target
     * like "7:5" with max_func=4 was parsed, stored, then never looked
     * at -- and the run printed "SUCCESS: Target distribution met!" over
     * a network containing no degree-7 nodes at all. Failing here makes
     * the request's impossibility visible instead of silently dropping
     * part of it. Mirrors the fail-fast guard in
     * PythonTopologyGenerator._validate_targets_reachable. */
    for (int i = max_func + 1; i <= max_possible_degree; ++i) {
        if (target_counts[i] > 0) {
            fprintf(stderr,
                    "Error: degree_distribution %d:%d is unreachable. Sculpting "
                    "enforces max_func=%d, so no node can finish with degree %d.\n",
                    i, target_counts[i], max_func, i);
            free(target_counts);
            freeGraph(base_graph);
            return 1;
        }
    }

    if (target_edge_count != -1) { // --- NEW ---
        printf("INFO: Target total edge count 'e' is set to %d.\n", target_edge_count);
    }

    // --- NEW: Validation for e:N ---
    if (target_edge_count != -1) {
        long long target_degree_sum = (long long)target_edge_count * 2;
        if (explicit_degree_sum > target_degree_sum) {
            fprintf(stderr, "Error: The sum of degrees from explicit targets (%lld) is already greater than the target total degree sum (%lld from e:%d).\n",
                    explicit_degree_sum, target_degree_sum, target_edge_count);
            free(target_counts);
            freeGraph(base_graph);
            return 1;
        }
        if (has_unconstrained) {
            // 'd:*' is the old "unconstrained" syntax, which is ambiguous with e:N
            fprintf(stderr, "Warning: Using e:%d (target edge count) overrides any 'd:*' (unconstrained) targets.\n", target_edge_count);
        }
    }
    // --- END NEW VALIDATION ---
    
    if (node_sum > (long long)base_graph->V) {
        fprintf(stderr, "Error: Sum of specified node counts in distribution (%lld) exceeds total nodes in %s lattice (%d).\n", node_sum, lattice_type, base_graph->V);
        free(target_counts);
        freeGraph(base_graph);
        return 1;
    }
    
    // --- MODIFIED: Only run handshake check if e:N is NOT set
    if (target_edge_count == -1 && !has_unconstrained && (explicit_degree_sum % 2 != 0)) {
        fprintf(stderr, "Error (Handshake Lemma): Sum of specified degrees (%lld) is odd.\n", explicit_degree_sum);
        free(target_counts);
        freeGraph(base_graph);
        return 1;
    }
    
    /* srand moved above the lattice build so MIX can draw its sites. */
    long long success_count = 0;

    printf("\n--- Initial State ---\n");
    print_distribution("Initial Lattice", base_graph, target_counts, -1, 0, max_func, extensive_logging);

    for (long long trial = 0; trial < max_trials; ++trial) {
        if (success_count >= max_saves) break;

        printf("\n--- Starting Trial %lld / %d (Found %lld so far) ---\n", trial, max_trials, success_count);

        // --- MODIFIED: Pass target_edge_count ---
        Graph* result_graph = run_single_trial(base_graph, max_func, target_counts, target_edge_count, trial, extensive_logging, dims_str, lattice_type);

        if (result_graph != NULL) {
            success_count++;
            printf("[Trial %lld | SUCCESS] Target distribution met!\n", trial);
            // MODIFIED: Pass dims_str (argv[1]) instead of N
            save_graph_to_file(result_graph, dims_str, trial);
            freeGraph(result_graph);
        } else {
            printf("[Trial %lld | FAILED] Could not find a valid network.\n", trial);
        }
    }
    
    printf("\n\nSIMULATION FINISHED.\n");
    if(success_count >= max_saves) {
        printf("Termination condition met: Found %lld networks (target was %d).\n", success_count, max_saves);
    } else {
        printf("All trials completed: Found %lld networks (target was %d).\n", success_count, max_saves);
    }
    
    freeGraph(base_graph);
    free(target_counts);
    return 0;
}