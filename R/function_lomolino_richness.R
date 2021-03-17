# Sugiro uma alternativa pra melhorar as estimativas de riqueza que s?o utilizadas como resposta no seu modelo. 
# Vc poderia utilizar uma fun??o param?trica com teto m?ximo pra aproximar a curva de rarefa??o.
# Imagina que vc gerou a curva de rarefa??o para um determinado s?tio (baseada em amostras) usando o pacote vegan, 
# e vc chamou esse objeto de rarefaction. A? vc cria um data.frame de duas colunas a partir desse output:

# n?o lembro se esses seriam os exatos nomes, mas s? pra ilustrar a ideia

#vegan_data <- data.frame(ln_richness = log(rarefaction$richness), samples = rarefaction$samples)

#Dps vc pode rodar o seguinte c?digo:

gen_lomolino_inits <- function(data) {

  data$richness <- exp(data$ln_richness)

  mean_inits <- getInitial(richness ~ SSlomolino(sample, Asym,

                                                 xmid, slope),

                           data = data)

  inits <- vector(mode = "list", length = nc)

  for (i in seq_along(inits)) {

    inits[[i]] <- lapply(mean_inits, function(x) {

      rnorm(1, mean = x, sd = abs(log10(abs(x)))) %>%

        abs %>%

        log

    })

    names(inits[[i]]) <- paste0("ln", tolower(names(inits[[i]])))

  }

  inits

}

##

run_lomolino_model <- function(lomolino_data) {

  formula <- brms::bf(ln_richness ~ lnasym -

                        log(1 + (exp(lnslope) ^ log(exp(lnxmid) / sample))),

                      lnasym + lnslope + lnxmid ~ 1,

                      nl = TRUE)

  my_inits <- gen_lomolino_inits(lomolino_data)

  priors <- 
            brms::prior(normal(0, 1), nlpar = "lnasym") +

            brms::prior(normal(0, 1), nlpar = "lnslope") +

            brms::prior(normal(0, 1), nlpar = "lnxmid")

  
  
  brms::brm(formula,

            data = lomolino_data, family = gaussian(),

            prior = priors, inits = my_inits,

            chains = nc,
            
            iter = ni,
            
            warmup = nb,
            
            thin=nt
            
  )

}


#model <- run_lomolino_model(vegan_data)

#Desse modelo voc? pode extrair a estimativa m?dia de lnasym, que representa a riqueza m?xima estimada daquele s?tio.

#asymp <- exp(fixefs(model)["lnasym_Intercept", "Estimate"])

# A fun??o Lomolino que codifiquei ? baseada no paper do Dengler et al. 2009 (veja a tabela 2)