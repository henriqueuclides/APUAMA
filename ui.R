# Define UI
ui <- fluidPage(theme = shinytheme("flatly"),
                useShinyjs(),
                tags$head(
                  tags$style(HTML(
                    "html {
                     position: relative;
                     min-height: 100%;
                   }
                   body {
                     height: 100%;
                     margin-bottom: 160px;
                     overflow-x: hidden;
                   }
                   .footer {
                     position: absolute;
                     left: 0px;
                     right: 0px;
                     bottom: 0;
                     width: 100%;
                     background-color: #2c3e50;
                     color: #fff;
                   }"))),
                tags$head(
                  tags$style(HTML("
                     .shiny-output-error-errorClass {
                      color: red;
                    }
                  "))
                ),
                
                navbarPage(
                  "APUAMA",
                  tabPanel("Project", icon = icon("atom"), align = "center",
                           fluidRow(style = "color: #fff; background-color: #1a242f",
                                    column(12,
                                           htmlOutput("Welcome"))
                           ),
                           htmlOutput("Intro"),
                           tags$img(src="workflow.jpg", height = "90%", width = "90%"),
                  ),
                  tabPanel("Tool", icon = icon("toolbox"),
                           sidebarPanel(
                             tags$h3("Input:"),
                             HTML(paste("An example can be downloaded", downloadLink("DownloadExamples", "here."),"<br/>&nbsp")),
                             fileInput("File1", "Choose ODS file for species"),
                             fileInput("File2", "Choose ODS file for reactions"),
                             actionButton("Showinputs","Show inputs"),
                             actionButton("Calculate","Calculate"),
                             hidden(tags$div(
                               id = "OBJoutput",
                               tags$h3("Output:"),
                               downloadButton("DownloadData","Download all data"),
                             ))
                           ),
                           mainPanel(
                             tags$div(
                               id = "OBJLets",
                               align="center",
                               tags$img(src="lets-begin.png", height = "90%", width = "90%")
                             ),
                             hidden(tags$div(
                               id = "OBJspecies",
                               align="center",
                               tags$h3("Species:"),
                               htmlOutput("Species", style = "font-weight: bold")
                             )),
                             hidden(tags$div(
                               id = "OBJreactions",
                               align="center",
                               tags$h3("Reactions:"),
                               htmlOutput("Reactions", style = "font-weight: bold")
                             )),
                             hidden(tags$div(
                               id = "OBJerrorinputs",
                               align="center",
                               tags$h3("Error:"),
                               htmlOutput("ErrorInputs", style = "color: red; font-weight: bold")
                             )),
                             hidden(tags$div(
                               id = "OBJerrorcalc",
                               tags$h3("Error:"),
                               htmlOutput("ErrorCalc", style = "color: red; font-weight: bold")
                             )),
                             hidden(tags$div(
                               id = "OBJinputs",
                               tabsetPanel(
                                 tabPanel("Species",
                                          div(
                                            style = "display:-webkit-flex; display:-ms-flexbox; display:flex;",
                                            uiOutput("SpeciesID"),
                                            hidden(tags$div(
                                              id = "OBJfileRov",
                                              style = "display:-webkit-flex; display:-ms-flexbox; display:flex;", div(style = "width: 30px;"),
                                              fileInput("File3", "Choose TXT file for rovibrational levels"), div(style = "width: 20px;"),
                                              htmlOutput("RovLevels", style = "font-weight: bold; font-size: 110%")
                                            )),
                                          ),
                                          hidden(tags$div(
                                            id = "OBJRov",
                                            tabsetPanel(
                                              tabPanel("Specie info",
                                                       htmlOutput("Inspecie"),
                                                       #div(dataTableOutput("Inspecie"), style = "font-size:80%")
                                              ),
                                              tabPanel(HTML("Table (<span>&#957</span>, J)"),
                                                       tags$h4(HTML("Table with the vibrational and rotational energy levels (<span>&#957</span>, J) in <i>kcal mol<sup>-1</sup></i>")),
                                                       div(dataTableOutput("InRvl", width = "100%", height = "auto"), style = "font-size:80%")
                                              ),
                                              tabPanel("Spectroscopic Constants",
                                                       htmlOutput("SpecConst"),
                                              ),
                                              tabPanel(HTML("Add (<span>&#957</span>, J) energy"),
                                                       div(
                                                         style = "display:flex;",
                                                         numericInput("levelv", label = HTML("<span>&#957</span>"), value = 0, min = 0, max = 49, width = "70px"),
                                                         div(style = "width: 10px;"),
                                                         numericInput("levelJ", label = "J", value = 0, min = 0, max = 49, width = "70px"),div(style = "width: 30px;"),
                                                         htmlOutput("NewEnergy")
                                                       ),
                                                       actionButton("EnterRov","Add energy")
                                              )
                                            )
                                          )),
                                          tags$div(
                                            id = "OBJtableSpecies",
                                            htmlOutput("Inspecies2"),
                                            #div(dataTableOutput("Inspecies"), style = "font-size:80%")
                                          ),
                                 ),
                                 tabPanel("Reactions",
                                          div(dataTableOutput("Inreactions"), style = "font-size:80%")
                                 )
                               )
                             )),
                             hidden(tags$div(
                               id = "OBJcalculate",
                               div(
                                 style = "display:-webkit-flex; display:-ms-flexbox; display:flex;",
                                 uiOutput("ReactionID"), div(style = "width: 30px;"), uiOutput("ReactionPrint")
                               ),
                               tabsetPanel(
                                 tabPanel("Reaction rate", align="center",
                                          uiOutput("Arrh"),
                                          checkboxInput('checkNa', HTML("divide by N<sub>A</sub>"), FALSE),
                                          plotOutput("PlotRate", width = "100%", height = "480")
                                 ),
                                 tabPanel("MEP",
                                          plotOutput("PlotMEP", width = "100%", height = "480")
                                 ),
                                 navbarMenu("Thermodynamic Properties",
                                            tabPanel("Enthalpy",
                                                     uiOutput("EqEao"),
                                                     plotOutput("PlotEao", width = "100%", height = "480")
                                            ),
                                            tabPanel("Entropy",
                                                     uiOutput("EqSto"),
                                                     plotOutput("PlotSto", width = "100%", height = "480")
                                            ),
                                            tabPanel("Heat Capacity",
                                                     uiOutput("EqCpo"),
                                                     plotOutput("PlotCpo", width = "100%", height = "480")
                                            )
                                 )
                               ),
                             )),
                           ),
                  ),
                  tabPanel("Instructions", icon = icon("align-justify"),
                           htmlOutput("Instructions"),
                           tabsetPanel(
                             tabPanel("Species Input", style = "color: ghostwhite; background-color: #1a242f",
                                      htmlOutput("InstrucSpecies")
                             ),
                             tabPanel("Reactions Input", style = "color: ghostwhite; background-color: #1a242f",
                                      htmlOutput("InstrucReactions")
                             ),
                             tabPanel("Rovibrational Levels Input", style = "color: ghostwhite; background-color: #1a242f",
                                      htmlOutput("InstrucRov")
                             )
                           ),
                  ),
                ),
                tags$footer(fluidRow(
                  column(1,),
                  column(5,
                         htmlOutput("Footer1")
                  ),
                  column(5,
                         htmlOutput("Footer2")
                  ),
                  column(1,)
                ), class = "footer"),
                setBackgroundColor("ghostwhite")
)
