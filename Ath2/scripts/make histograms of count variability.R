#f is the outcome of filter_counts_based_on_treatment()
ft <- f %>% t() %>% as.data.frame() %>% gather ()
ggplot(ft, aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x') +
  scale_y_log10()

f_no_0 <- f[vapply(f,
                                          function(z) length(unique(z)) > 1,
                                          logical(1L))] %>% data.matrix(., rownames.force = T)

median_all_genes <- median(f_no_0)
median_per_gene <- apply(f_no_0, 2, median)
f_no_0_t <- f_no_0 %>% t() %>% as.data.frame()
f_f <- f_no_0_t %>% filter(median_per_gene > median_all_genes) %>% t() %>% as.data.frame()

f_ft <- f_f %>% t() %>% as.data.frame() %>% gather()
head(f_ft)
ggplot(f_ft, aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x') +
  scale_y_log10() 
