#ifndef DIXON_WIN_GUI_WRAPPER_H
#define DIXON_WIN_GUI_WRAPPER_H

#include <stddef.h>

typedef enum {
    DIXON_GUI_MODE_SOLVE = 0,
    DIXON_GUI_MODE_RESULTANT = 1,
    DIXON_GUI_MODE_COMPLEXITY = 2,
    DIXON_GUI_MODE_IDEAL = 3
} dixon_gui_mode_t;

typedef struct {
    dixon_gui_mode_t mode;
    const char *dixon_path;
    const char *polynomials;
    const char *eliminate_vars;
    const char *field;
    const char *ideal_generators;
    const char *omega;
    int solve_verbose;
} dixon_gui_request_t;

typedef struct {
    int exit_code;
    char *resolved_dixon_path;
    char *command_line;
    char *stdout_text;
    char *stderr_text;
    char *combined_output;
} dixon_gui_response_t;

void dixon_gui_response_init(dixon_gui_response_t *response);
void dixon_gui_response_clear(dixon_gui_response_t *response);

int dixon_gui_find_default_cli(char *buffer, size_t buffer_size);
int dixon_gui_run_request(const dixon_gui_request_t *request,
                          dixon_gui_response_t *response,
                          char *error_message,
                          size_t error_message_size);

#endif
