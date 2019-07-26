# Code of example 7
#
# Works under Windows
#
#
#

# Set the RNG seed
rng_seed <- 314
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  arg <- suppressWarnings(as.numeric(args[1]))
  if (is.na(arg)) {
    stop(
      "Please supply a numerical value for the RNG seed. \n",
      "Actual value: ", args[1]
    )
  }
  rng_seed <- arg
  if (rng_seed < 1) {
    stop("Please supply an RNG seed with a positive non-zero value")
  }
}
if (length(args) > 1) {
  stop(
    "Please supply only 1 argument for the RNG seed. \n",
    "Number of arguments given: ", length(args) - 1
  )
}

library(pirouette)
suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))
library(beautier)

root_folder <- getwd()
example_no <- 7
example_folder <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
setwd(example_folder)
set.seed(rng_seed)
testit::assert(is_beast2_installed())
phylogeny <- create_yule_tree(n_taxa = 6, crown_age = 10)

alignment_params <- create_alignment_params(
  root_sequence = create_blocked_dna(length = 1000),
  mutation_rate = 0.1,
  rng_seed = rng_seed
)

experiment <- create_gen_experiment()
experiments <- list(experiment)

# Set the RNG seed
for (i in seq_along(experiments)) {
  experiments[[i]]$beast2_options$rng_seed <- rng_seed
}

# Testing
if (1 == 2) {
  for (i in seq_along(experiments)) {
    experiments[[i]]$inference_model$mcmc <- create_mcmc(chain_length = 20000, store_every = 1000)
  }
}

pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments,
  twinning_params = create_twinning_params(
    rng_seed_twin_tree = rng_seed,
    rng_seed_twin_alignment = rng_seed
  )
)

################################################################################
# Settings to run on Peregrine cluster
################################################################################
pir_params$alignment_params$fasta_filename <- file.path(example_folder, "true.fasta")
for (i in seq_along(pir_params$experiments)) {
  pir_params$experiments[[i]]$beast2_options$input_filename <- file.path(example_folder, "beast2_input.xml")
  pir_params$experiments[[i]]$beast2_options$output_log_filename <- file.path(example_folder, "beast2_output.log")
  pir_params$experiments[[i]]$beast2_options$output_trees_filenames <- file.path(example_folder, "beast2_output.trees")
  pir_params$experiments[[i]]$beast2_options$output_state_filename <- file.path(example_folder, "beast2_output.xml.state")
  pir_params$experiments[[i]]$beast2_options$beast2_working_dir <- example_folder
  pir_params$experiments[[i]]$errors_filename <- file.path(example_folder, "error.csv")
  pir_params$experiments[[i]]$beast2_options$overwrite <- TRUE
}
pir_params$evidence_filename <- file.path(example_folder, "evidence_true.csv")
if (!is_one_na(pir_params$twinning_params)) {
  pir_params$twinning_params$twin_tree_filename <- file.path(example_folder, "twin.tree")
  pir_params$twinning_params$twin_alignment_filename <- file.path(example_folder, "twin.fasta")
  pir_params$twinning_params$twin_evidence_filename <- file.path(example_folder, "evidence_twin.csv")
}
rm_pir_param_files(pir_params)
################################################################################

errors <- pir_run(
  phylogeny,
  pir_params = pir_params
)

if (1 == 2) {
  errors <- utils::read.csv(
    file = file.path(example_folder, "errors.csv")
  )
  check_pir_out(errors)

  pir_plot(errors)
}

utils::write.csv(
  x = errors,
  file = file.path(example_folder, "errors.csv"),
  row.names = FALSE
)


pir_plot(errors) +
  ggsave(file.path(example_folder, "errors.png"))

pir_to_pics(
  phylogeny = phylogeny,
  pir_params = pir_params,
  folder = example_folder
)

pir_to_tables(
  pir_params = pir_params,
  folder = example_folder
)
