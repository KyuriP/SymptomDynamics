
### plot just for the intro ppt:
sde_out <- euler_stochastic2(
  deterministic_rate = dif_eq,
  stochastic_rate = sto_eq,
  initial_condition = init,
  parameters1 = parms1H,
  deltaT = deltaT,
  timelength = 5000,
  D1 = 0.01,
  D2 = D_stoeq1,
  shock = TRUE,
  parameters2 = parms2H, 
  t_shock = t_shock,
  duration = shock_duration,
  seed = 12#6789
)

tot <- sde_out |>
  mutate(totalsymptom = rowSums(pick(S_anh:S_sui))) 

t_shock <- 1000
shock_duration <- 1000

shock_period <- data.frame(time = 0:timelength) |>
  mutate(shock = ifelse(time >= t_shock & time <= t_shock + shock_duration, TRUE, FALSE))


hysplot <- tot |>
  ggplot(aes(x = t, y = totalsymptom)) +
  geom_line(col = "red4")+
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])*ncol(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  scale_y_continuous(breaks = c(0, 8), labels = c("Healthy", "Depressed"))+
  labs(x = "time", y= "", title="Total Symptom Activation") +
  theme_classic(base_size = 40)+
  theme(axis.title.x = element_text(vjust = 9, hjust = 1.05),
        plot.margin = unit(c(0.5, 2, 0, 0),  "inches")) 


ggsave("hysplot.png", plot = hysplot, width = 80, height =30, units = "cm", dpi = 300)
ggsave("hysplot_high.png", plot = hysplot, width = 80, height =30, units = "cm", dpi = 300)


### resilience
x <- seq(-2, 2, by = 0.01)

anh_ref = x^2*(22116*x^3 - 26745*x^2- 10760*x + 14340) / 12000
anh_low = x^2*(1552*x^3 - 1915*x^2- 830*x + 1195) / 1000
anh_high = x^2*(12804*x^3 - 15255*x^2- 5780*x + 7170) / 6000



anh_landscapes <- data.frame(x, anh_ref, anh_low, anh_high) |>
  pivot_longer(!x, values_to = "val", names_to = "resilience") |>
  mutate(resilience = factor(resilience, level = c("anh_low", "anh_ref", "anh_high")))


anh_land <- anh_landscapes |>
  ggplot() +
  geom_line(aes(x = x, y = val, col = resilience), linewidth = 0.6) +
  scale_color_manual(values= c("royalblue","darkorange", "coral3"), labels = c("high", "base", "low"))+
  theme_classic()+
  ylim(-0.2, 0.3) + 
  labs(color = "Resilience", y ="", x = "Anhedonia") +
  scale_x_continuous(breaks=c(0,1), limits =c(-0.2, 1.2))+
  #guides(fill = guide_legend(byrow = TRUE)) +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        plot.title=element_text(size = 20), 
        axis.title.x=element_text( size = 20), 
        axis.title.y=element_text(hjust = 1, vjust = 1, size = 20, angle = 360),
        axis.text.x=element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face ="bold"),
        legend.spacing.y = unit(10, "pt"),
        legend.key.width = unit(3, "line"),
        legend.position='top')

ggsave("anh_resilience.png", plot = anh_land, width = 30, height =20, units = "cm", dpi = 300)
