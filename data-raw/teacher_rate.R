## code to prepare `teacher_rate` dataset

dat <- readr::read_csv("../zhzh-text-data/prof1000.stem.gender.csv")
rate_prof <- na.omit(dat)
rate_prof <- rate_prof %>%
  select(rating, grade, comments)

docs_vocab <- prep_docs(na.omit(dat), "comments")
docsdf <- as_tibble(docs_vocab$documents, .name_repair = "unique")

rate_wdocs <- rate_prof %>%
  select(-comments) %>%
  bind_cols(docsdf)

# Really rely on the bag-of-words assumption in topic models
rate_syn <- synthpop::syn(rate_wdocs, minnumlevels = 5, seed = 333)

rate_syn2 <- rate_syn$syn %>%
  mutate(...76 = as.numeric(as.character(...76)),
         ...77 = as.numeric(as.character(...76)),
         id = row_number()) %>%
  tidyr::pivot_longer(starts_with("..."), names_to = "position", values_to = "index") %>%
  select(id, rating, grade, position, index) %>%
  rowwise() %>%
  mutate(index = if_else(near(index, 0), NA_real_, index)) %>%
  ungroup()

vocab_df <- tibble(word = docs_vocab$vocab) %>%
  mutate(index = row_number())

teacher_rate <- rate_syn2 %>%
  left_join(vocab_df, by = "index") %>%
  select(-index) %>%
  tidyr::pivot_wider(id_cols = id:position, names_from = "position", values_from = "word") %>%
  tidyr::unite(doc, starts_with("..."), sep = " ", na.rm = TRUE)

usethis::use_data(teacher_rate, overwrite = TRUE, compress = "bzip2")
