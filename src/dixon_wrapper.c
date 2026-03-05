// dixon_wrapper.c - Fully isolated wrapper module containing all FLINT-related code
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Cross-platform DLL/shared library export definition  __declspec(dllexport)
#ifdef DLL_EXPORT
	#define EXPORT
#else
	#define EXPORT
#endif

// Include all FLINT headers - safe to load after DLL path is configured
#include <flint/fmpz_factor.h>
#include "dixon_flint.h"
#include "dixon_interface_flint.h"
#include "fq_mvpoly.h"
#include "fq_unified_interface.h"
#include "unified_mpoly_resultant.h"
#include "dixon_with_ideal_reduction.h"
#include "polynomial_system_solver.h"

// Internal prime power check function
static int check_prime_power_internal(const fmpz_t n, fmpz_t prime, unsigned long *power) {
    if (fmpz_cmp_ui(n, 1) <= 0) return 0;
    
    // Check if it is directly a prime
    if (fmpz_is_probabprime(n)) {
        fmpz_set(prime, n);
        *power = 1;
        return 1;
    }
    
    // Use FLINT's factorization functionality
    fmpz_factor_t factors;
    fmpz_factor_init(factors);
    fmpz_factor(factors, n);
    
    // Check if there is only one prime factor (i.e. a prime power)
    if (factors->num == 1) {
        fmpz_set(prime, factors->p + 0);  // First (and only) prime factor
        *power = factors->exp[0];         // Corresponding exponent
        fmpz_factor_clear(factors);
        return 1;
    }
    
    fmpz_factor_clear(factors);
    return 0;
}

// Internal function to parse field size
static int parse_field_size_internal(const char *field_str, fmpz_t prime, unsigned long *power) {
    if (!field_str || strlen(field_str) == 0) {
        return 0;
    }
    
    // Check if the string contains '^' (power notation)
    const char *caret = strchr(field_str, '^');
    if (caret) {
        // Format: p^k
        char *prime_str = malloc(caret - field_str + 1);
        strncpy(prime_str, field_str, caret - field_str);
        prime_str[caret - field_str] = '\0';
        
        // Parse the prime part
        fmpz_t p;
        fmpz_init(p);
        int success = fmpz_set_str(p, prime_str, 10);
        if (success != 0) {
            fmpz_clear(p);
            free(prime_str);
            return 0;
        }
        
        // Parse the exponent part
        char *endptr;
        unsigned long k = strtoul(caret + 1, &endptr, 10);
        if (*endptr != '\0' || k == 0) {
            fmpz_clear(p);
            free(prime_str);
            return 0;
        }
        
        // Check if p is prime
        if (!fmpz_is_probabprime(p)) {
            fmpz_clear(p);
            free(prime_str);
            return 0;
        }
        
        fmpz_set(prime, p);
        *power = k;
        fmpz_clear(p);
        free(prime_str);
        return 1;
    } else {
        // Format: plain number (may be p or p^k)
        fmpz_t field_size;
        fmpz_init(field_size);
        int success = fmpz_set_str(field_size, field_str, 10);
        if (success != 0) {
            fmpz_clear(field_size);
            return 0;
        }
        
        // Check if it is a prime power
        int result = check_prime_power_internal(field_size, prime, power);
        fmpz_clear(field_size);
        return result;
    }
}

// Public interface: validate field size
EXPORT int validate_field_size(const char *field_str, char *error_msg, int error_msg_size) {
    fmpz_t prime;
    unsigned long power;
    fmpz_init(prime);
    
    int result = parse_field_size_internal(field_str, prime, &power);
    
    if (!result) {
        if (error_msg && error_msg_size > 0) {
            snprintf(error_msg, error_msg_size, "Failed to parse '%s' as a valid field size", field_str);
        }
    } else {
        // Check if the prime is too large
        if (!fmpz_fits_si(prime)) {
            if (error_msg && error_msg_size > 0) {
                snprintf(error_msg, error_msg_size, "Prime is too large (must fit in machine word)");
            }
            result = 0;
        }
    }
    
    fmpz_clear(prime);
    return result;
}

// Public interface: get field info
EXPORT char* get_field_info(const char *field_str) {
    fmpz_t prime;
    unsigned long power;
    fmpz_init(prime);
    
    if (!parse_field_size_internal(field_str, prime, &power)) {
        fmpz_clear(prime);
        return strdup("Invalid field");
    }
    
    char *info = malloc(256);
    if (!info) {
        fmpz_clear(prime);
        return NULL;
    }
    
    if (power == 1) {
        // Prime field
        char *prime_str = fmpz_get_str(NULL, 10, prime);
        snprintf(info, 256, "Field: F_%s (prime field)", prime_str);
        //free(prime_str);
    } else {
        // Extension field
        unsigned long prime_ui = fmpz_get_ui(prime);
        unsigned long field_size = 1;
        for (unsigned long i = 0; i < power; i++) {
            field_size *= prime_ui;
        }
        snprintf(info, 256, "Field: F_%lu^%lu (extension field, size %lu, generator 't')", 
                prime_ui, power, field_size);
    }
    
    fmpz_clear(prime);
    return info;
}

// Public interface: basic Dixon resultant computation
EXPORT char* dixon_compute_basic(const char *polys_str, const char *vars_str, const char *field_str) {
    fmpz_t prime;
    unsigned long power;
    fmpz_init(prime);
    
    if (!parse_field_size_internal(field_str, prime, &power)) {
        fmpz_clear(prime);
        return strdup("Error: Invalid field size");
    }
    
    // Initialize finite field context
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, prime, power, "t");
    
    // Call Dixon resultant function
    char *result = dixon_str(polys_str, vars_str, ctx);
    
    // Cleanup
    //fq_nmod_ctx_clear(ctx);
    //fmpz_clear(prime);
    
    //if (!result) {
        //return strdup("Error: Dixon resultant computation failed");
    //}
    
    return result;
}


// Simple solution formatting in computation thread - GCC COMPATIBLE VERSION WITH WINDOWS LINE ENDINGS
static char* format_solutions_simple(const polynomial_solutions_t *sols, const fq_nmod_ctx_t ctx) {
    if (!sols) {
        return strdup("Solution structure is null\r\n");
    }
    
    char *buffer = malloc(16384);  // 16KB buffer
    if (!buffer) {
        return strdup("Memory allocation failed\r\n");
    }
    
    strcpy(buffer, "\r\n=== POLYNOMIAL SYSTEM SOLUTIONS ===\r\n");
    
    if (!sols->is_valid) {
        strcat(buffer, "Solving failed");
        if (sols->error_message) {
            strcat(buffer, ": ");
            strncat(buffer, sols->error_message, 500);
        }
        strcat(buffer, "\r\n");
        return buffer;
    }
    
    if (sols->has_no_solutions == 1) {
        strcat(buffer, "System has no solutions over the finite field\r\n");
        return buffer;
    }
    
    if (sols->has_no_solutions == -1) {
        strcat(buffer, "Solving failed: polynomial system dimension greater than zero, please use Dixon resultant elimination\r\n");
        return buffer;
    }
    
    if (sols->num_variables == 0) {
        strcat(buffer, "No variables\r\n");
        return buffer;
    }
    
    if (sols->num_solution_sets == 0) {
        strcat(buffer, "No solutions found\r\n");
        return buffer;
    }
    
    char temp[512];
    sprintf(temp, "Found %ld complete solution set(s):\r\n", sols->num_solution_sets);
    strcat(buffer, temp);
    
    // Format solutions with simple error handling
    for (slong set = 0; set < sols->num_solution_sets && set < 5; set++) {
        sprintf(temp, "\r\nSolution set %ld:\r\n", set + 1);
        strcat(buffer, temp);
        
        for (slong var = 0; var < sols->num_variables && var < 10; var++) {
            if (!sols->variable_names || !sols->variable_names[var]) {
                continue;
            }
            
            sprintf(temp, "  %s = ", sols->variable_names[var]);
            strcat(buffer, temp);
            
            if (!sols->solutions_per_var) {
                strcat(buffer, "ERROR: solutions_per_var is null\r\n");
                continue;
            }
            
            slong num_sols = sols->solutions_per_var[set * sols->num_variables + var];
            
            if (num_sols == 0) {
                strcat(buffer, "no solution\r\n");
            } else if (num_sols == 1) {
                if (!sols->solution_sets || !sols->solution_sets[set] || 
                    !sols->solution_sets[set][var] || !sols->solution_sets[set][var][0]) {
                    strcat(buffer, "ERROR: Invalid solution pointer\r\n");
                } else {
                    // Simple approach: try conversion, use fallback if it fails
                    char *sol_str = NULL;
                    
                    // Try the pretty print conversion
                    sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set][var][0], ctx);
                    
                    // Check if conversion was successful
                    if (sol_str && strlen(sol_str) > 0 && strlen(sol_str) < 500) {
                        strncat(buffer, sol_str, 200);
                        strcat(buffer, "\r\n");
                    } else {
                        // Fallback: for prime fields, try to display as a simple integer
                        if (fq_nmod_ctx_degree(ctx) == 1) {
                            // For prime fields, extract the coefficient
                            fmpz_t coeff;
                            fmpz_init(coeff);
                            fq_nmod_get_fmpz(coeff, sols->solution_sets[set][var][0], ctx);
                            
                            if (fmpz_fits_si(coeff)) {
                                slong val = fmpz_get_si(coeff);
                                sprintf(temp, "%ld\r\n", val);
                                strcat(buffer, temp);
                            } else {
                                strcat(buffer, "[large integer value]\r\n");
                            }
                            fmpz_clear(coeff);
                        } else {
                            strcat(buffer, "[field extension element]\r\n");
                        }
                    }
                    
                    // Clean up sol_str if it was allocated
                    if (sol_str) {
                        ;//free(sol_str); // This causes a segfault
                    }
                }
            } else if (num_sols > 1 && num_sols <= 5) {
                strcat(buffer, "{");
                for (slong sol = 0; sol < num_sols; sol++) {
                    if (sol > 0) strcat(buffer, ", ");
                    
                    if (sols->solution_sets && sols->solution_sets[set] && 
                        sols->solution_sets[set][var] && sols->solution_sets[set][var][sol]) {
                        
                        // For multiple solutions, use simple fallback
                        if (fq_nmod_ctx_degree(ctx) == 1) {
                            fmpz_t coeff;
                            fmpz_init(coeff);
                            fq_nmod_get_fmpz(coeff, sols->solution_sets[set][var][sol], ctx);
                            
                            if (fmpz_fits_si(coeff)) {
                                slong val = fmpz_get_si(coeff);
                                sprintf(temp, "%ld", val);
                                strcat(buffer, temp);
                            } else {
                                strcat(buffer, "[large_val]");
                            }
                            fmpz_clear(coeff);
                        } else {
                            strcat(buffer, "[ext_field]");
                        }
                    } else {
                        strcat(buffer, "NULL");
                    }
                }
                strcat(buffer, "}\r\n");
            } else {
                sprintf(temp, "%ld solutions (too many to display)\r\n", num_sols);
                strcat(buffer, temp);
            }
        }
    }
    
    strcat(buffer, "\r\n=== Solution Complete ===\r\n");
    return buffer;
}

// Public interface: polynomial system solver
EXPORT char* dixon_compute_solver(const char *polys_str, const char *field_str) {
    fmpz_t prime;
    unsigned long power;
    fmpz_init(prime);
    
    if (!parse_field_size_internal(field_str, prime, &power)) {
        fmpz_clear(prime);
        return strdup("Error: Invalid field size");
    }
    
    // Initialize finite field context
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, prime, power, "t");
    
    // Call polynomial system solver
    polynomial_solutions_t *solutions = solve_polynomial_system_string(polys_str, ctx);
    char *result = format_solutions_simple(solutions, ctx);
    // Cleanup
    //fq_nmod_ctx_clear(ctx);
    //fmpz_clear(prime);
    
    //if (!result) {
        //return strdup("Error: Polynomial system solving failed");
    //}
    
    return result;
}

// Public interface: Dixon with ideal reduction computation
EXPORT char* dixon_compute_ideal(const char *polys_str, const char *vars_str, 
                         const char *ideal_str, const char *field_str) {
    fmpz_t prime;
    unsigned long power;
    fmpz_init(prime);
    
    if (!parse_field_size_internal(field_str, prime, &power)) {
        fmpz_clear(prime);
        return strdup("Error: Invalid field size");
    }
    
    // Initialize finite field context
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, prime, power, "t");
    
    // Call Dixon with ideal reduction function
    char *result = dixon_with_ideal_reduction_str(polys_str, vars_str, ideal_str, ctx);
    
    // Cleanup
    fq_nmod_ctx_clear(ctx);
    fmpz_clear(prime);
    
    if (!result) {
        return strdup("Error: Dixon resultant with ideal reduction failed");
    }
    
    return result;
}
