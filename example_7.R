# Code of example 7
#
# Works under Windows
library(pirouette)
suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))
library(beautier)

root_folder <- getwd()
example_no <- 7
example_folder <- file.path(root_folder, paste0("example_", example_no))
dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
setwd(example_folder)
set.seed(314)
testit::assert(is_beast2_installed())
phylogeny <- create_yule_tree(n_taxa = 6, crown_age = 10)

alignment_params <- create_alignment_params(
  root_sequence = create_blocked_dna(length = 1000),
  mutation_rate = 0.1
)

experiment <- create_gen_experiment()
experiments <- list(experiment)

# Testing
if (1 == 2) {
  for (i in seq_along(experiments)) {
    experiments[[i]]$inference_model$mcmc <- create_mcmc(chain_length = 20000, store_every = 1000)
  }
}

pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments,
  twinning_params = create_twinning_params()
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

testit::assert(pir_params$experiments[[1]]$inference_model$mcmc$store_every != -1)
esses <- tracerer::calc_esses(
  traces = tracerer::parse_beast_log(pir_params$experiments[[1]]$beast2_options$output_log_filename),
  sample_interval = pir_params$experiments[[1]]$inference_model$mcmc$store_every
)

df_esses <- data.frame(parameter = colnames(esses), ESS = as.character(esses))

sink(file.path(example_folder, "esses.latex"))
xtable::print.xtable(
  xtable::xtable(df_esses, caption = "ESSes of example 1", label = "tab:esses_example_1", digits = 0),
  include.rownames = FALSE
)
sink()

print("#######################################################################")
print("Appendix")
print("#######################################################################")
pir_to_pics(
  phylogeny = phylogeny,
  pir_params = pir_params,
  folder = example_folder
)
