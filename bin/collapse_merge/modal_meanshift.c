#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_ITER 100

typedef struct {
    int key;  // original data point
    int val;  // count or assigned mode
} IntMap;

// Parse compressed input: <int>:<count>,<int>:<count>,...
IntMap *parse_input(const char *input_str, int *n) {
    char *input = strdup(input_str);
    char *token = strtok(input, ",");
    int capacity = 128;
    IntMap *data = malloc(sizeof(IntMap) * capacity);
    int count = 0;

    while (token) {
        if (count >= capacity) {
            capacity *= 2;
            data = realloc(data, sizeof(IntMap) * capacity);
            if (!data) {
                fprintf(stderr, "Memory allocation failed\n");
                exit(1);
            }
        }
        char *colon = strchr(token, ':');
        if (!colon) {
            fprintf(stderr, "Invalid input format: %s\n", token);
            exit(1);
        }
        *colon = '\0';
        data[count].key = atoi(token);
        data[count].val = atoi(colon + 1);
        count++;
        token = strtok(NULL, ",");
    }

    free(input);
    *n = count;
    return data;
}

// Find weighted mode in the window
int find_mode_in_window(IntMap *data, int n, double center, double bandwidth) {
    int mode = -1;
    int max_weight = -1;

    for (int i = 0; i < n; i++) {
        if (fabs(data[i].key - center) <= bandwidth) {
            if (data[i].val > max_weight || (data[i].val == max_weight && data[i].key < mode)) {
                max_weight = data[i].val;
                mode = data[i].key;
            }
        }
    }

    return mode;
}

void mean_shift_modes(IntMap *data, int n, double bandwidth, double epsilon) {
    for (int i = 0; i < n; i++) {
        double x = data[i].key;
        for (int iter = 0; iter < MAX_ITER; iter++) {
            int new_mode = find_mode_in_window(data, n, x, bandwidth);
            if (fabs(new_mode - x) < epsilon) break;
            x = new_mode;
        }
        printf("%d:%d", data[i].key, (int)x);
        if (i != n - 1) printf(",");
    }
    printf("\n");
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <input_file> <bandwidth> <epsilon>\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    double bandwidth = atof(argv[2]);
    double epsilon = atof(argv[3]);

    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Failed to open input file");
        return 1;
    }

    char *line = NULL;
    size_t len = 0;

    while (getline(&line, &len, fp) != -1) {
        // Remove newline characters
        line[strcspn(line, "\r\n")] = '\0';

        // Split line into ID and data string
        char *tab = strchr(line, '\t');
        if (!tab) {
            fprintf(stderr, "Invalid line format: %s\n", line);
            continue;
        }
        *tab = '\0';
        char *id = line;
        char *data_str = tab + 1;

        int n;
        IntMap *data = parse_input(data_str, &n);

        printf("%s\t", id);
        mean_shift_modes(data, n, bandwidth, epsilon);

        free(data);
    }

    free(line);
    fclose(fp);
    return 0;
}

