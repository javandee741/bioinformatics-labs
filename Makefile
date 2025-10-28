# === Paths ===
PYTHON=python
SCRIPTS=scripts
RESULTS=results
FILES=files
INFRA=infra

# === Targets ===

# Анализ аминокислот (задание недели 5)
amino-stats:
	$(PYTHON) $(SCRIPTS)/amino_counter.py \
		-i $(FILES)/protein.fasta \
		-o $(RESULTS)/amino_stats.tsv \
		--plots-dir $(RESULTS)

# Сборка Docker-образа
docker-build:
	docker build -t bioinfo:week5 -f $(INFRA)/Dockerfile.bio .

# Запуск контейнера
docker-run:
	docker run --rm \
		-v $(PWD)/files:/app/files \
		-v $(PWD)/results:/app/results \
		bioinfo:week5

# Запуск пайплайна Snakemake (для будущих заданий)
snakemake-run:
	snakemake --cores 2

# Очистка результатов
clean:
	rm -f $(RESULTS)/*.tsv $(RESULTS)/*.png

#Snakemake
snakemake-run:
	snakemake --cores 2 --printshellcmds

snakemake-force:
	snakemake --cores 2 --forceall --printshellcmds

snakemake-rerun:
	snakemake --cores 2 -R amino_stats --printshellcmds

clean-results:
	rm -f results/*.tsv results/*.png

