setwd("~/git/DA-analysis")
library(ArchR)
library(tidyverse)
library(magrittr)
library(argparse)

parser = ArgumentParser(prog = 'inner-create-ArchR-project.R')
grid = read.delim("sh/grids/preprocessing/create-ArchR-project.txt")
for (param_name in colnames(grid))
	parser$add_argument(paste0('--', param_name),
						type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

if (!dir.exists(args$archR_dir)){
	dir.create(args$archR_dir, recursive=T)
}
setwd(args$archR_dir)

addArchRGenome(args$ref_genome)
addArchRThreads(threads = 8)

input_files = list.files(args$fragments_dir, full.names=T)
reformatFragmentFiles(input_files)	

names(input_files) = gsub('_fragments.*', '', basename(input_files))

ArrowFiles = createArrowFiles(
	inputFiles = input_files,
	sampleNames = names(input_files),
	minTSS = 0, #Dont set this too high because you can always increase later
	minFrags = 0, 
	addTileMat = TRUE,
	addGeneScoreMat = TRUE
)

proj = ArchRProject(
	ArrowFiles = ArrowFiles, 
	outputDirectory = "final_proj",
	copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

saveArchRProject(ArchRProj = proj, outputDirectory = "final_proj", load = FALSE)