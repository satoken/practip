/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.6
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "PRactIP"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "PRactIP"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "0.0.1"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  int threads_arg;	/**< @brief The number of threads for IP solver (default='1').  */
  char * threads_orig;	/**< @brief The number of threads for IP solver original value given at command line.  */
  const char *threads_help; /**< @brief The number of threads for IP solver help description.  */
  char * train_arg;	/**< @brief Train the parameters from given data.  */
  char * train_orig;	/**< @brief Train the parameters from given data original value given at command line.  */
  const char *train_help; /**< @brief Train the parameters from given data help description.  */
  char * predict_arg;	/**< @brief Predict interactions.  */
  char * predict_orig;	/**< @brief Predict interactions original value given at command line.  */
  const char *predict_help; /**< @brief Predict interactions help description.  */
  int cross_validation_arg;	/**< @brief Perform the n-fold cross validation (default='0').  */
  char * cross_validation_orig;	/**< @brief Perform the n-fold cross validation original value given at command line.  */
  const char *cross_validation_help; /**< @brief Perform the n-fold cross validation help description.  */
  float eta_arg;	/**< @brief Initial step width for the subgradient optimization (default='0.5').  */
  char * eta_orig;	/**< @brief Initial step width for the subgradient optimization original value given at command line.  */
  const char *eta_help; /**< @brief Initial step width for the subgradient optimization help description.  */
  float pos_w_arg;	/**< @brief The weight for positive interactions (default='4').  */
  char * pos_w_orig;	/**< @brief The weight for positive interactions original value given at command line.  */
  const char *pos_w_help; /**< @brief The weight for positive interactions help description.  */
  float neg_w_arg;	/**< @brief The weight for negative interactions (default='1').  */
  char * neg_w_orig;	/**< @brief The weight for negative interactions original value given at command line.  */
  const char *neg_w_help; /**< @brief The weight for negative interactions help description.  */
  float reg_w_arg;	/**< @brief The weight for the L1 regularization term (default='0.125').  */
  char * reg_w_orig;	/**< @brief The weight for the L1 regularization term original value given at command line.  */
  const char *reg_w_help; /**< @brief The weight for the L1 regularization term help description.  */
  float semi_w_arg;	/**< @brief The weight for the graph regularization term for semi-supervised learning (default='1.0').  */
  char * semi_w_orig;	/**< @brief The weight for the graph regularization term for semi-supervised learning original value given at command line.  */
  const char *semi_w_help; /**< @brief The weight for the graph regularization term for semi-supervised learning help description.  */
  int d_max_arg;	/**< @brief The maximim number of iterations of the supervised learning (default='25').  */
  char * d_max_orig;	/**< @brief The maximim number of iterations of the supervised learning original value given at command line.  */
  const char *d_max_help; /**< @brief The maximim number of iterations of the supervised learning help description.  */
  int g_max_arg;	/**< @brief The maximum number of iterations of the semi-supervised learning (default='5').  */
  char * g_max_orig;	/**< @brief The maximum number of iterations of the semi-supervised learning original value given at command line.  */
  const char *g_max_help; /**< @brief The maximum number of iterations of the semi-supervised learning help description.  */
  int aa_int_max_arg;	/**< @brief The maximum number of interations of each amino acid (default='3').  */
  char * aa_int_max_orig;	/**< @brief The maximum number of interations of each amino acid original value given at command line.  */
  const char *aa_int_max_help; /**< @brief The maximum number of interations of each amino acid help description.  */
  int rna_int_max_arg;	/**< @brief The maximum number of interations of each nucleotide (default='4').  */
  char * rna_int_max_orig;	/**< @brief The maximum number of interations of each nucleotide original value given at command line.  */
  const char *rna_int_max_help; /**< @brief The maximum number of interations of each nucleotide help description.  */
  float exceeding_penalty_arg;	/**< @brief The penalty for exceeding the limit of the number of interactions for each residue/base (default='0.5').  */
  char * exceeding_penalty_orig;	/**< @brief The penalty for exceeding the limit of the number of interactions for each residue/base original value given at command line.  */
  const char *exceeding_penalty_help; /**< @brief The penalty for exceeding the limit of the number of interactions for each residue/base help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int threads_given ;	/**< @brief Whether threads was given.  */
  unsigned int train_given ;	/**< @brief Whether train was given.  */
  unsigned int predict_given ;	/**< @brief Whether predict was given.  */
  unsigned int cross_validation_given ;	/**< @brief Whether cross-validation was given.  */
  unsigned int eta_given ;	/**< @brief Whether eta was given.  */
  unsigned int pos_w_given ;	/**< @brief Whether pos-w was given.  */
  unsigned int neg_w_given ;	/**< @brief Whether neg-w was given.  */
  unsigned int reg_w_given ;	/**< @brief Whether reg-w was given.  */
  unsigned int semi_w_given ;	/**< @brief Whether semi-w was given.  */
  unsigned int d_max_given ;	/**< @brief Whether d-max was given.  */
  unsigned int g_max_given ;	/**< @brief Whether g-max was given.  */
  unsigned int aa_int_max_given ;	/**< @brief Whether aa-int-max was given.  */
  unsigned int rna_int_max_given ;	/**< @brief Whether rna-int-max was given.  */
  unsigned int exceeding_penalty_given ;	/**< @brief Whether exceeding-penalty was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
