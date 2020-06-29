# Paquetes
library(ggplot2)
library(dplyr)
library(reshape2)
library(magrittr)
library(ggiraph)

# Para fuentes
library(extrafont)
#font_import()
loadfonts(device = "win")

# Datos: Global Cancer Observatory (World Health Organization)
stats <- read.csv2("marimekko_gco_mortalidad.csv", encoding = "UTF-8")

casos_otros <- c(stats[nrow(stats), "Hombres"] - sum(stats[1:(nrow(stats) - 1), "Hombres"]),
           stats[nrow(stats), "Mujeres"] - sum(stats[1:(nrow(stats) - 1), "Mujeres"])) %>% as.numeric

# Se quita la categorías del total y se añade la de "Otros"
stats <- rbind(stats[-nrow(stats), ], c("Otras localizaciones", casos_otros)) %>% 
  melt(id.vars = "Loc")

names(stats) <- c("localizacion", "Sexo", "casos")
stats$casos <- as.numeric(stats$casos)
stats$localizacion <- ordered(factor(stats$localizacion, levels = unique(stats$localizacion)))

# Referencia: https://dqn.website/post/interactive-mekko-charts-in-r/
stats %<>% group_by(localizacion) %>%
  mutate(
    share = casos / sum(casos),
    tot_group = sum(casos)
  ) %>% ungroup()

stats %<>%
  group_by(Sexo) %>% 
  arrange(desc(localizacion)) %>%
  mutate(
    ymax = cumsum(tot_group) / sum(tot_group), 
    ymin = (ymax - (tot_group/sum(tot_group)))
  ) %>% ungroup() %>% 
  group_by(localizacion) %>% 
  arrange(desc(Sexo)) %>%
  mutate(xmax = cumsum(share), xmin = xmax - share) %>%
  ungroup() %>% 
  arrange(localizacion)

# Etiquetas para el gráfico
labels <- stats %>% 
  mutate(y = ymax - 0.01, yRange = (ymax - ymin)* 100) %>%
  mutate(label = ifelse(casos != 0, paste0(localizacion, " - ", Sexo, " (", format(casos, trim = TRUE, big.mark = "."), " casos)"), NA)) %>% 
  select(label, xmin, y, yRange) %>% 
  ungroup()

print(labels)

ggplot(stats) + 
  geom_rect(aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, fill = Sexo),
            colour = "white") +
  geom_text(data = labels,
            aes(x = xmin + 0.008, y = y, label = label), 
            hjust = 0, vjust = 1, colour = "white", size = 5, family = "Perpetua") +
  scale_x_continuous(position = "top", expand = c(0.01, 0.01), 
                     labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.02)) +
  scale_fill_manual(values = c("turquoise4", "steelblue3")) +
  ggtitle("Mortalidad por cáncer en el mundo, 2018",
          "Distribución de defunciones por sexo y localización anatómica. \nFuente: Global Cancer Observatory, Organización Mundial de la Salud.") +
  theme(axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 14, color = "black", family = "Perpetua", face = "bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 24, hjust = .5, family = "Perpetua"),
        plot.subtitle = element_text(face = "italic", size = 18, hjust = .5, family = "Perpetua"),
        legend.position = "none",
        panel.background = element_blank())

ggsave("marimekko_gco_mortalidad.png", width = 10, height = 10, dpi = 600)

