plan <- drake_plan(
  ## load data
  speed = read_csv("data/video_swimming_speeds.csv"),
  respiro = read_csv("data/data_respirometry.csv"),
  speedmax = read_csv("data/data_fulton_2007.csv"),
  moorea = read_csv("data/data_uvc_moo_2016_fsp.csv"),
  ## wrangle data
  speedmax_complete = wrangle_speedmax(speedmax),
  ## run models
  mod_mr = run_mod_mr(respiro),
  mod_speed = run_mod_speed(speed),
  mod_speedmax = run_mod_speedmax(speedmax_complete),
  ## get R2
  mod_mr_R2 = bayes_R2(mod_mr),
  mod_speed_R2 = bayes_R2(mod_speed),
  mod_speedmax_R2 = bayes_R2(mod_speedmax),
  ## get tables
  table_mr = get_table_mr(mod_mr),
  table_speed = get_table_speed(mod_speed),
  #table_speedmax = get_table_speedmax(mod_speedmax),
  ## FMR calculation
  ### 1) reference dataframe
  reference = create_reference(moorea),
  ### 2) get fmr, fsa, fas
  field_summary = get_fmr(reference, mod_speed, mod_speedmax, mod_mr),
  ### 3) get metabolic scaling 
  slopes = get_slopes(field_summary, mod_mr),
  ## figures
  plot_mr = make_plot_mr(mod_mr, respiro),
  plot_speed = make_plot_speed(reference, mod_speed, mod_speedmax, speed),
  plot_fsa = make_plot_fsa(field_summary),
  plot_slopes = make_plot_slopes(slopes),
  plot_combined = make_plot_combined(plot_fsa, plot_slopes),
  #plot_community_mr = make_plot_community_mr(field_summary, moorea),
  plot_community_ab = make_plot_community_ab(moorea),
  figure5 = make_fig5(moorea, field_summary),
  plot_annex = make_plot_sizecor(),
  fig1 = load_fig1(),
  ## text
  main_text_doc = rmarkdown::render(knitr_in("text/main_text.Rmd"), 
                                    output_format = "word_document", 
                                    output_dir = "./output/text/",
                                    output_file = "main.docx")
)
