-- create view for CNV handling  
CREATE MATERIALIZED VIEW cnv_view AS
SELECT d.id, s.code AS subject, m.chr->>'value' AS chr, (m.start->>'value')::int AS start, (m.stop->>'value')::int AS stop,
        m.cytoband_start->>'value' AS cytoband_start, m.cytoband_stop->>'value' AS cytoband_stop,
        (m.stop->>'value')::int - (m.start->>'value')::int AS base_length,
        int8range((m.start->>'value')::int, (m.stop->>'value')::int, '[]') AS base_interval,
        (m.amplification->>'value')::float AS amplification, (m.deletion->>'value')::float AS deletion, 
        (m.pval->>'value')::float AS pvalue, 
        ARRAY(SELECT jsonb_array_elements_text(m.gene_name->'values')) AS genes, 
        ARRAY(SELECT jsonb_array_elements_text(m.mirna->'values')) AS mirnas
FROM data d
LEFT JOIN subject s ON s.id = d.parent_subject
CROSS JOIN LATERAL jsonb_to_record(d.metadata) AS m(chr jsonb, start jsonb, stop jsonb, cytoband_start jsonb, cytoband_stop jsonb, 
                                                    amplification jsonb, deletion jsonb, pval jsonb, gene_name jsonb, mirna jsonb)
WHERE d.type = 8;

CREATE VIEW cgh_genomic_profile_view AS
SELECT d.id, s.code, d.metadata->'type'->>'value' AS profile_type, d.metadata->'sca_type'->>'value' AS sca_type
FROM data d
LEFT JOIN subject s ON s.id = d.parent_subject
WHERE d.type = 18;

CREATE INDEX cnv_chr ON cnv_view USING btree(chr);
CREATE INDEX cnv_start ON cnv_view USING btree(start);
CREATE INDEX cnv_stop ON cnv_view USING btree(stop);
CREATE INDEX cnv_base_length ON cnv_view USING btree(base_length);
CREATE INDEX cnv_amplification ON cnv_view USING btree(amplification);
CREATE INDEX cnv_deletion ON cnv_view USING btree(deletion);
CREATE INDEX cnv_base_interval_idx ON cnv_view USING gist(base_interval);


-- count aberration calls per band interval (amplifications)
SELECT count(*), chr, cytoband_start, cytoband_stop from cnv_view WHERE amplification > 0 GROUP BY (chr, cytoband_start, cytoband_stop) ORDER BY count(*) DESC;

-- count  aberration calls per band interval (deletions)
SELECT count(*), chr, cytoband_start, cytoband_stop from cnv_view WHERE deletion < 0 GROUP BY (chr, cytoband_start, cytoband_stop) ORDER BY count(*) DESC;


-- group genes by occurences in the CNVS (amplifications)
SELECT array_agg(res.genes), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(genes) AS genes FROM cnv_view WHERE deletion < 0 GROUP BY genes) AS res GROUP BY count ORDER BY count DESC;

-- group genes by occurences in the CNVS (amplifications)
SELECT array_agg(res.genes), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(genes) AS genes FROM cnv_view WHERE amplification > 0 GROUP BY genes) AS res GROUP BY count ORDER BY count DESC;

-- group mirnas by occurences in the CNVS (amplifications)
SELECT array_agg(res.mirnas), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(mirnas) AS mirnas FROM cnv_view WHERE deletion < 0 GROUP BY mirnas) AS res GROUP BY count ORDER BY count DESC;

-- group mirnas by occurences in the CNVS (amplifications)
SELECT array_agg(res.mirnas), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(mirnas) AS mirnas FROM cnv_view WHERE amplification > 0 GROUP BY mirnas) AS res GROUP BY count ORDER BY count DESC;

-- group genes by occurences in the CNVS (amplifications)
SELECT string_agg(res.genes, ','), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(genes) AS genes FROM cnv_view WHERE deletion < 0 GROUP BY genes) AS res GROUP BY count ORDER BY count DESC;

-- group genes by occurences in the CNVS (amplifications)
SELECT string_agg(res.genes, ','), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(genes) AS genes FROM cnv_view WHERE amplification > 0 GROUP BY genes) AS res GROUP BY count ORDER BY count DESC;

-- group mirnas by occurences in the CNVS (amplifications)
SELECT string_agg(res.mirnas, ','), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(mirnas) AS mirnas FROM cnv_view WHERE deletion < 0 GROUP BY mirnas) AS res GROUP BY count ORDER BY count DESC;

-- group mirnas by occurences in the CNVS (amplifications)
SELECT string_agg(res.mirnas, ','), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(mirnas) AS mirnas FROM cnv_view WHERE amplification > 0 GROUP BY mirnas) AS res GROUP BY count ORDER BY count DESC;

/* COPY commands */
\copy (SELECT count(*), chr, cytoband_start, cytoband_stop from cnv_view WHERE deletion < 0 GROUP BY (chr, cytoband_start, cytoband_stop) ORDER BY count(*) DESC) TO 'cytoband_occurrences_del.csv' DELIMITER ';' CSV HEADER

\copy (SELECT count(*), chr, cytoband_start, cytoband_stop from cnv_view WHERE amplification > 0 GROUP BY (chr, cytoband_start, cytoband_stop) ORDER BY count(*) DESC) TO 'cytoband_occurrences_amp.csv' DELIMITER ';' CSV HEADER

PREPARE search_cnv_by_gene(int, jsonb) AS
WITH s AS (
    SELECT id, code, sex, personal_info FROM subject
), pd AS (
    SELECT id, given_name, surname, birth_date FROM personal_details
) 
SELECT DISTINCT d.id, s.code, s.sex, pd.given_name, pd.surname, pd.birth_date, d.metadata FROM data d 
LEFT JOIN s ON s.id = d.parent_subject 
LEFT JOIN pd ON pd.id = s.personal_info 
WHERE d.type = $1 AND ((d.metadata @> $2));

\copy (SELECT string_agg(res.genes, ','), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(genes) AS genes FROM cnv_view WHERE deletion < 0 GROUP BY genes) AS res GROUP BY count ORDER BY count DESC) TO 'gene_occurrences_del.csv' DELIMITER ';' CSV HEADER

\copy (SELECT string_agg(res.genes, ','), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(genes) AS genes FROM cnv_view WHERE amplification > 0 GROUP BY genes) AS res GROUP BY count ORDER BY count DESC) TO 'gene_occurrences_amp.csv' DELIMITER ';' CSV HEADER

\copy (SELECT string_agg(res.mirnas, ','), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(mirnas) AS mirnas FROM cnv_view WHERE deletion < 0 GROUP BY mirnas) AS res GROUP BY count ORDER BY count DESC) TO 'mirna_occurrences_del.csv' DELIMITER ';' CSV HEADER

\copy (SELECT string_agg(res.mirnas, ','), res.count FROM (SELECT count(*) AS count, jsonb_array_elements_text(mirnas) AS mirnas FROM cnv_view WHERE amplification > 0 GROUP BY mirnas) AS res GROUP BY count ORDER BY count DESC) TO 'mirna_occurrences_amp.csv' DELIMITER ';' CSV HEADER

-- PREPARED STATEMENT for searching genes in CNV
PREPARE search_cnv_by_gene(int, jsonb) AS
WITH s AS (
    SELECT id, code, sex, personal_info FROM subject
), pd AS (
    SELECT id, given_name, surname, birth_date FROM personal_details
)   
SELECT DISTINCT d.id, s.code, s.sex, pd.given_name, pd.surname, pd.birth_date, d.metadata FROM data d
LEFT JOIN s ON s.id = d.parent_subject
LEFT JOIN pd ON pd.id = s.personal_info
WHERE d.type = $1 AND ((d.metadata @> $2));


-- PL/pgSQL queries
CREATE OR REPLACE FUNCTION test_cnv() RETURNS integer AS $$
DECLARE 
    gene_curs cursor FOR SELECT jsonb_array_elements_text(genes) AS genes FROM cnv_view GROUP BY genes LIMIT 5000;
    gene_text text;
    gene_json json;
    gene_jsonb jsonb;
    cnv_type integer := 8;
    rows integer;
BEGIN
    FOR gene IN gene_curs 
    LOOP
        gene_text := '{"gene_name":{"values":["' || gene.genes || '"]}}';
        gene_json := to_json(gene_text);
        gene_jsonb := gene_json::jsonb;
        RAISE NOTICE '%', gene_jsonb;
        SELECT count(*) INTO rows FROM data d
            WHERE d.type = cnv_type AND ((d.metadata @> gene_jsonb));
        RAISE NOTICE 'Found % rows', rows;
    END LOOP;
    RETURN 0;
END;
$$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION count_genes_in_cnv() RETURNS integer AS $$
DECLARE
    cnvs cursor FOR SELECT id, ARRAY(SELECT jsonb_array_elements_text(genes)) AS genes FROM cnv_view;
    cnv_len integer;
    total_len integer := 0;
BEGIN
    FOR cnv in cnvs LOOP
        SELECT INTO cnv_len cardinality(cnv.genes);
        RAISE NOTICE 'genes array has for CNV % length %', cnv.id, cnv_len;
        total_len = total_len + cnv_len;
        RAISE NOTICE 'total length: %', total_len;
    END LOOP;
    RETURN total_len;
END 
$$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION compute_gene_occurrences_in_cnv() RETURNS TABLE(gene text, count integer) AS $$
DECLARE
    cnvs cursor FOR SELECT id, ARRAY(SELECT jsonb_array_elements_text(genes)) AS genes FROM cnv_view;
    genes_arr text[] := '{}';
    total_len integer := 0;
BEGIN
    FOR cnv in cnvs LOOP
        genes_arr = array_cat(genes_arr, cnv.genes);
    END LOOP;
    RETURN QUERY gene, count(*) FROM ( SELECT unnest(genes_arr) ) AS gene GROUP BY gene ORDER BY count(*) DESC;
END 
$$ LANGUAGE plpgsql;


CREATE VIEW genes_view AS
SELECT id, ARRAY (SELECT jsonb_array_elements_text(genes)) AS genes, cardinality(ARRAY (SELECT jsonb_array_elements_text(genes))) AS genes_count
FROM cnv_view;

CREATE OR REPLACE VIEW genes_view AS
SELECT id, ARRAY(SELECT jsonb_array_elements_text(genes)) AS genes, cardinality(ARRAY (SELECT jsonb_array_elements_text(genes))) AS genes_count
FROM cnv_view;

-- GET single gene occurrences
WITH sub1 AS (
    SELECT genes FROM genes_view
)
SELECT res AS gene, count(*) FROM (SELECT unnest(genes) AS res FROM sub1) AS gene GROUP BY gene ORDER BY count(*) DESC;

WITH sub1 AS (
    SELECT genes FROM genes_view
),
sub2 as (
    SELECT res AS gene, count(*) FROM (SELECT unnest(genes) AS res FROM sub1) AS gene GROUP BY gene ORDER BY count(*) DESC
)
SELECT count(*) FROM sub2;

CREATE OR REPLACE VIEW mirnas_view AS
SELECT id, ARRAY(SELECT jsonb_array_elements_text(mirnas)) AS mirnas, cardinality(ARRAY (SELECT jsonb_array_elements_text(mirnas))) AS mirnas_count
FROM cnv_view;

WITH sub1 AS (
    SELECT mirnas FROM mirnas_view
),
sub2 as (
    SELECT res AS mirnas, count(*) FROM (SELECT unnest(mirnas) AS res FROM sub1) AS mirnas GROUP BY mirnas ORDER BY count(*) DESC
)
SELECT count(*) FROM sub2;

-- Microdeletions export
\copy (SELECT * FROM cnv_view WHERE base_length < 1000000 AND deletion < 0 ORDER BY id) TO 'microdeletions.csv' DELIMITER ';' CSV HEADER

\copy (SELECT id, subject, chr, start, stop, cytoband_start, cytoband_stop, deletion, pvalue, array_to_string(genes, ',') AS genes, array_to_string(mirnas, ',') AS mirnas FROM cnv_view WHERE base_length < 1000000 AND deletion < 0 ORDER BY id) TO 'microdeletions.csv' DELIMITER ';' CSV HEADER


-- Regular Human Genome 19 variants table
-- NOTE: rename 'end' (SQL keyword) header to 'stop' 
CREATE TABLE grch37_hg19_variants (
    id serial PRIMARY KEY,
    variant_accession text UNIQUE NOT NULL,
    chr text NOT NULL,
    start integer NOT NULL,
    stop integer NOT NULL,
    base_interval int8range,
    base_length integer,
    variant_type text NOT NULL,
    variant_subtype text NOT NULL,
    reference text NOT NULL,
    pubmed_id integer NOT NULL,
    method text NOT NULL,
    platform text,
    merged_variants text,
    supporting_variants text NOT NULL,
    merged_or_sample text NOT NULL,
    frequency float,
    sample_size integer NOT NULL,
    observed_gains integer,
    observed_losses integer,
    cohort_description text,
    genes text,
    samples text,
    supporting_variants_arr text[],
    genes_arr text[],
    samples_arr text[]
);

CREATE UNIQUE INDEX grch37_variant_accession_idx ON grch37_hg19_variants USING btree(variant_accession);
CREATE INDEX grch37_chr ON grch37_hg19_variants USING btree(chr);
CREATE INDEX grch37_start ON grch37_hg19_variants USING btree(start);
CREATE INDEX grch37_stop ON grch37_hg19_variants USING btree(stop);
CREATE INDEX grch37_base_interval_idx ON grch37_hg19_variants USING gist(base_interval);
CREATE INDEX grch37_variant_subtype ON grch37_hg19_variants USING btree(variant_subtype);

-- TRIGGER
CREATE OR REPLACE FUNCTION handle_hg19_variants_on_insert() RETURNS TRIGGER AS $handle_hg19_variants_on_insert$
BEGIN
    NEW.chr := concat('chr', NEW.chr);
    NEW.base_interval := int8range(NEW.start, NEW.stop, '[]');
    NEW.supporting_variants_arr := string_to_array(NEW.supporting_variants, ',', '');
    NEW.genes_arr := string_to_array(NEW.genes, ',', '');
    NEW.samples_arr := string_to_array(NEW.samples, ',', '');
    RETURN NEW;
END;
$handle_hg19_variants_on_insert$ LANGUAGE plpgsql;

CREATE TRIGGER grch37_hg19_variants_trigger 
    BEFORE INSERT ON grch37_hg19_variants
    FOR EACH ROW
    EXECUTE PROCEDURE handle_hg19_variants_on_insert(); 

\copy grch37_hg19_variants (variant_accession, chr, start, stop, variant_type, variant_subtype, reference, pubmed_id, method, platform, merged_variants, supporting_variants, merged_or_sample, frequency, sample_size, observed_gains, observed_losses, cohort_description, genes, samples) from '~/Projects/unige/GRCh37_hg19_variants_2016-05-15.csv' with delimiter ';' csv header

-- convert TSV to semicolon-CSV with
sed $'s/\t/;/g' GRCh37_hg19_variants_2016-05-15.txt >  GRCh37_hg19_variants_2016-05-15.csv

SELECT c.id, c.chr, c.start, c.stop, v.variant_accession
FROM cnv_view c
    LEFT JOIN grch37_hg19_variants v ON c.chr = v.chr AND c.base_interval && v.base_interval;

\copy (SELECT c.id AS call_id, c.chr AS call_chr, c.start AS call_start, c.stop AS call_stop, v.variant_accession AS variant_accession, v.chr AS var_chr, v.start AS var_start, v.stop AS var_stop FROM cnv_view c LEFT JOIN grch37_hg19_variants v ON c.chr = v.chr AND c.base_interval && v.base_interval WHERE c.deletion < 0 AND c.base_length <= 1000000 AND v.variant_subtype = 'loss' ORDER BY id) TO 'microdels-calls-1MB.csv' DELIMITER ';' CSV HEADER

SELECT c.id, c.subject, c.chr, c.start, c.stop, c.base_length, c.cytoband_start, c.cytoband_stop, c.deletion, c.pvalue, string_agg(v.variant_accession, ',') AS variants 
FROM cnv_view c 
    LEFT JOIN grch37_hg19_variants v ON c.chr = v.chr AND c.base_interval && v.base_interval 
WHERE c.deletion < 0 AND c.base_length <= 1000000 AND v.variant_subtype = 'loss' 
GROUP BY c.id, c.subject, c.chr, c.start, c.stop, c.base_length, c.cytoband_start, c.cytoband_stop, c.deletion, c.pvalue
ORDER BY c.id;

-- ANNOTATE the CNVs with the variants contained in 
\copy (SELECT c.id, c.subject, c.chr, c.start, c.stop, c.base_length, c.cytoband_start, c.cytoband_stop, c.deletion, c.pvalue, string_agg(v.variant_accession, ',') AS variants FROM cnv_view c LEFT JOIN grch37_hg19_variants v ON c.chr = v.chr AND c.base_interval && v.base_interval WHERE c.deletion < 0 AND c.base_length <= 5000000 AND v.variant_subtype = 'loss' GROUP BY c.id, c.subject, c.chr, c.start, c.stop, c.base_length, c.cytoband_start, c.cytoband_stop, c.deletion, c.pvalue ORDER BY c.id) to 'microdels-calls-5MB.csv' delimiter ';' csv header

-- Query to find overlapping variants
SELECT c.id, c.subject, c.chr, c.start, c.stop, c.base_length, c.cytoband_start, c.cytoband_stop, c.deletion, c.pvalue, string_agg(v.variant_accession, ',') AS variants 
FROM cnv_view c
    LEFT JOIN grch37_hg19_variants v ON c.chr = v.chr AND c.base_interval && v.base_interval 
WHERE EXISTS (
    SELECT 1 FROM cnv_view c2
    WHERE c.chr = c2.chr AND c.base_interval && c2.base_interval AND c.base_length < 5000000 AND c.deletion < 0
)
AND v.variant_subtype = 'loss'
GROUP BY c.id, c.subject, c.chr, c.start, c.stop, c.base_length, c.cytoband_start, c.cytoband_stop, c.deletion, c.pvalue
ORDER by chr, start ASC;
