# Define server functions
server <- function(input, output) {
  therespecies <<- FALSE
  therereactions <<- FALSE
  
  # Welcome
  output$Welcome <- renderText({
    "<h2><b>Welcome kinetics researchers!"
  })
  
  # Intro
  output$Intro <- renderText({
    "<p style=\"text-align:justify; font-size:110%;\"><br/>APUAMA is a user friendly R app designed to determine the reaction rate
      and the thermodynamic properties of uni/bimolecular reactions. With data from electronic structure calculations, the code,
      which is based on the transition state theory, determines the rate constant with tunneling corrections, such as Wigner, Eckart
      and small curvature, and can also include the rovibrational level of diatomic molecules. The results are presented in
      Arrhenius-Kooij form for the reaction rate, and the thermodynamic properties are written in polynomial form. APUAMA means “fast”
      in Tupi-Guarani language, so the code calculates the reaction rate in a simple and intuitive graphical interface, quickly
      and conveniently."
  })
  
  # set download link to examples
  output$DownloadExamples <- downloadHandler(
    filename = function(){
      paste0("examples.zip")},
    content = function(file){
      #go to a temp dir to avoid permission issues
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files <- NULL;
      
      df <- data.frame(name = c("","1-CO","","","2-CO","","","C","","CO2","","","","","OCOC","","","","",""),
                       natoms = c("",2,"","",2,"","",1,"",3,"","","","",4,"","","","",""),
                       nfreqs = c("",1,"","",1,"","",0,"",4,"","","","",6,"","","","",""),
                       symmetry = c("",1,"","",1,"","",1,"",2,"","","","",1,"","","","",""),
                       type = c("","linear","","","linear","","","atom","","linear","","","","","nonlinear","","","","",""),
                       energy = c("",-113.19037050,"","",-113.19037050,"","",-37.73550530,"",-188.38959510,"","","","",-226.07478620,"","","","",""),
                       atoms = c("","C","O","","C","O","","C","","C","O","O","","","C","O","C","O","",""),
                       x = c("",0.,0.,"",0.,0.,"",0.,"",0.,0.,0.,"","",0.,0.,1.3805067213,2.5241841525,"",""),
                       y = c("",0.,0.,"",0.,0.,"",0.,"",0.,0.,0.,"","",0.,0.,0.,-0.0000247879 ,"",""),
                       z = c("",0.,1.13185722,"",0.,1.13185722,"",0.,"",0.,1.16312786,-1.16312786,"","",0.,1.29542567,1.2173214011,0.8487068711,"",""),
                       frequencies = c("",2160.1425,"","",2160.1425,"","","","",669.9411,670.1048,1349.5657,2389.0185,"",-1485.7424,374.6883,438.9152,821.9531,1305.5782,1924.3383)
      )
      write_ods(df, "species.ods", row_names = FALSE, col_names = TRUE)
      files <- c("species.ods",files)
      
      df <- data.frame(reactant1 = c("","1-CO"),
                       reactant2 = c("","2-CO"),
                       ts = c("","OCOC"),
                       product1 = c("","C"),
                       product2 = c("","CO2"),
                       psi = c("",0)
      )
      write_ods(df, "reaction.ods", row_names = FALSE, col_names = TRUE)
      files <- c("reaction.ods",files)
      
      write("1-CO
 5.2577533758163115
 10.747286764903201
 16.01044859511485
 14.687050365930753
 5.805396776586387
 0.8889869546207309
 335.8122320350773
 1.1389824523677687"
            ,"1-co-rovlevels.txt",sep = ' ')
      files <- c("1-co-rovlevels.txt",files)
      
      write("2-CO
 5.2577533758163115
 10.747286764903201
 16.01044859511485
 14.687050365930753
 5.805396776586387
 0.8889869546207309
 335.8122320350773
 1.1389824523677687"
            ,"2-co-rovlevels.txt",sep = ' ')
      files <- c("2-co-rovlevels.txt",files)
      
      
      #create the zip file
      zip(file,files)
    }
  )
  
  # Instructions
  output$Instructions <- renderText({
    paste0("<p style=\"text-align:justify;\">Although the usability with APUAMA is quite simple, there are some details that must be taken care of. The code is based
             on TST, so each reaction must have a transition state, either for unimolecular or bimolecular reactions, and can be of the type:
             <br><b>A", icon("long-arrow-alt-right"), "B; &nbsp &nbsp A", icon("long-arrow-alt-right"), "B+C; &nbsp &nbsp
             A+B", icon("long-arrow-alt-right"), "C; &nbsp &nbsp A+B", icon("long-arrow-alt-right"), "C+D</b><br>
             where <b>A, B, C</b> and <b>D</b> can be atoms or molecules.<br>
             There are only two required inputs, one for <b>species</b> and one for <b>reactions</b>, and an optional input for the
             <b>rovibrational levels</b> calculation of diatomic molecules. <br>After the calculations, a .zip file will be available for download with all calculated data, such as the reaction rate as \"rateTSNAME.ods\",
             the minimum energy path  as \"mepTSNAME.ods\", the thermodynamic properties (enthalpy, entropy and heat capacity) and the polynomial
             fitting with seven and nine NASA coefficients for each specie will be written down as \"thermopropSPECIENAME.ods\".<br>")
  })
  
  # Instructions for species input
  output$InstrucSpecies <- renderText({
    df <- data.frame(name = c("","C","","H2","","","H2CO","","","","","","","TS","","","","",""),
                     natoms = c("",1,"",2,"","",4,"","","","","","",4,"","","","",""),
                     nfreqs = c("",0,"",1,"","",6,"","","","","","",6,"","","","",""),
                     symmetry = c("",1,"",2,"","",2,"","","","","","",1,"","","","",""),
                     type = c("","atom","","linear","","","nonlinear","","","","","","","nonlinear","","","","",""),
                     energy = c("",-37.7355053,"",-1.163732,"","", -114.341304,"","","","","","",-114.209903,"","","","",""),
                     atoms = c("","C","","H","H","","C","O","H","H","","","","C","O","H","H","",""),
                     x = c("",0.,"",0.,0.,"",0.132148,-0.037178,1.010045,-0.578527,"","","",0.297907,0.150231,0.257392,-0.546612,"",""),
                     y = c("",0.,"",0.,0.,"",-0.212138,0.062703,0.152644,-0.847963,"","","",-0.314778,-1.360162,0.750940,0.924000,"",""),
                     z = c("",0.,"",0.,0.744531,"",0.009596,1.165250,-0.563472,-0.558507,"","","",0.554056,1.031658,0.323088,1.356282,"",""),
                     frequencies = c("","","",4414.2230,"","",1202.3073,1270.3180,1539.1572,1827.2856,2868.8925,2918.1350,"",-1865.7806,789.6168,929.7337,1331.0581,1931.8362,3194.2344)
    )
    paste0("<div style=\"margin:10px;\">This entry must have a tabular format in the .ods extension. It will be necessary to include all the species
      that will be used, as well as the TSs.
             The following header is required:

             <br><br><b>name natoms nfreqs symmetry type energy atoms x y z frequencies</b><br><br>

             <p style=\"text-align:justify;\">Followed by the specie name, the number of atoms, the number of vibrational frequencies, an integer number for the external symmetry,
             the geometry type of the specie <font color=#18b598>(atom, linear or nonlinear)</font>, the calculated energy <font color=#18b598>(in Hartree)</font>,
             each defined atom in the molecule with its optimized geometry in cartesian coordinates <font color=#18b598>(in <span>&#8491</span>)</font>,
             and for last, the vibrational frequencies <font color=#18b598>(in cm<sup>-1</sup>)</font>.
             <br>Like the example:",
           htmlTable(df,caption = "&nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       "),
           "<p style=\"text-align:justify;\"><b>Notes: All species names must be different and each of then cannot exceed 50 characters.
             APUAMA supports calculations with atoms from H to Xe, and the name of each atom must be the same as in the periodic table (ex: H, He, Li, Be... etc).
             The order of species does not matter. The imaginary frequency for the TS should be included as its first frequency.
             If the specie does not have a frequency the respective cell must be empty. The rest of the table cells must be blank.")
  })
  
  # Instructions for reactions input
  output$InstrucReactions <- renderText({
    df <- data.frame(reactant1 = c("","H2CO","H2CO","H2"),
                     reactant2 = c("","ZERO","ZERO","CO"),
                     ts = c("","TS1","TS2","TS3"),
                     product1 = c("","H2OC","H2","H2O"),
                     product2 = c("","ZERO","CO","C"),
                     psi = c("",0,0,0)
    )
    paste0("<div style=\"margin:10px;\">This entry must have a tabular format in the .ods extension. It will be necessary to include all the reactions per line.
             The following header is required:

             <br><br><b>reactant1 reactant2 ts product1 product2 psi</b><br><br>

             <p style=\"text-align:justify;\">Followed by the species names for the first and second reactants, the TS, and the first and second products,
             the <b>psi</b> parameter is a correction factor for the vibrational frequencies, is a number close to 1, and has a default value of
             <font color=#18b598>0.9705646</font> if psi was set to <font color=#18b598>0</font>.<br>
             Like the example:",
           htmlTable(df,caption = "&nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                       "),
           "<p style=\"text-align:justify;\"><b>Notes: The order of the reactions does not matter. For reactions without one of the reactants or products, the word ZERO should be included,
             in this case the order matters, so the second reactant and/or the second product must be ZERO.")
  })
  
  # Instructions for RVL input
  output$InstrucRov <- renderText({
    "<div style=\"margin:10px;\"><p style=\"text-align:justify;\">
      This entry is a simple text file in the .txt extension. Is a hidden feature only available for diatomic molecules, to access it, all the others data entries (for the species and reactions) must be
      included first, then, access is given by pressing the button \"Show Inputs\" and choosing the respective diatomic molecule for data entry.<br>
      For the parameters, the name of the molecule should be included first (the name must be equal as in the species entry), the following data must have a blank space in the beggining of the line
      and each of them separated by a line break, such as the five parameters of the Generalized Rydberg potential function of fifth degree, the reference energy <font color=#18b598>(in kcal mol<sup>-1</sup>)</font>,
      the dissociation energy <font color=#18b598>(in kcal mol<sup>-1</sup>)</font> and the equilibrium distance <font color=#18b598>(in <span>&#8491</span>)</font>.<br>
             Like the example for CO:<br><br>

             CO<br>
             &nbsp 5.25775<br>
             &nbsp 10.74728<br>
             &nbsp 16.01044<br>
             &nbsp 14.68705<br>
             &nbsp 5.80539<br>
             &nbsp 0.88898<br>
             &nbsp 335.81223<br>
             &nbsp 1.13898<br>

             <p style=\"text-align:justify;\"> After entering the properly data, the rovibrational levels will be calculated and presented in a table with
      the respective vibrational and rotational energy levels and as spectroscopic constants, as well as the addition of their level that must be entered manually.
      <br><b>Notes: To reset the specie energy after the addition of a rovibrational level just include the level (0,0) and the energy will be restored.
      Each diatomic molecule that will be added to the rovibrational level should have its own data input."
  })
  
  # Our information footer1
  output$Footer1<-renderText({
    paste0("<div style=\"margin:10px;\"><h4>",icon("newspaper"),"&nbsp <b>APUAMA:</b> a software tool for reaction rate calculations</h4>
             Euclides, H.O., P. Barreto, P.R. J Mol Model 23, 176 (2017).<br>",
           tags$a(href="https://link.springer.com/article/10.1007/s00894-017-3337-5", "https://doi.org/10.1007/s00894-017-3337-5"),
           "<br>Henrique de Oliveira Euclides &nbsp", icon("envelope"), "<i>&nbsp henrique.euclides@inpe.br</i><br/>
             Patrícia Regina Pereira Barreto &nbsp", icon("envelope"), "<i>&nbsp patricia.barreto@inpe.br</i>")
  })
  
  # Our information footer2
  output$Footer2<-renderText({
    paste0("<div style=\"margin:10px;\">Laboratório Associado de Plasma–LAP, Instituto Nacional de Pesquisas Espaciais–INPE/MCT, CP515, São José dos Campos, SP CEP 12247-970, Brazil<br/><br/>",
           tags$img(src = "inpe-logo.png", width = "62px", height = "50px"),
           "&nbsp &nbsp &nbsp &nbsp &nbsp",
           tags$img(src = "capes-logo.png", width = "50px", height = "50px"))
  })
  
  # error message for entry
  output$ErrorInputs <- renderText({
    "Error: enter with species and reactions inputs"
  })
  
  # input for species
  output$Species <- reactive({
    therespecies <<- FALSE
    shinyjs::hide("OBJoutput")
    shinyjs::hide("OBJLets")
    file <- input$File1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "ods", "Error: please upload a .ods file for the species"), errorClass = "errorClass")
    
    species <<- read_ods(path =  file$datapath, col_names = TRUE)
    
    nspecies <<- sum(!is.na(species$natoms))
    namespecies <<- species$name[!is.na(species$name)]
    natmospecies <<- species$natoms[!is.na(species$natoms)]
    nfreqspecies <<- species$nfreqs[!is.na(species$nfreqs)]
    symmetryspecies <<- species$symmetry[!is.na(species$symmetry)]
    typespecies <<- species$type[!is.na(species$type)]
    energyspecies <<- species$energy[!is.na(species$energy)]
    atomspecies <<- species$atoms[!is.na(species$atoms)]
    xspecies <<- species$x[!is.na(species$x)]
    yspecies <<- species$y[!is.na(species$y)]
    zspecies <<- species$z[!is.na(species$z)]
    frequenciespecies <<- species$frequencies[!is.na(species$frequencies)]
    
    extractSpecies(nspecies,namespecies,natmospecies,nfreqspecies,symmetryspecies,typespecies,energyspecies,atomspecies,xspecies,yspecies,zspecies,frequenciespecies)
    
    # define which species to show
    output$SpeciesID <- renderUI({
      selectInput("speciesID", "Choose specie:", namespecies, width = '150')
    })
    
    # set the actual selected specie
    setSpecie <- function(){
      j = 1
      if(!is.null(input$speciesID)){
        for (i in 1:nspecies)
          if(input$speciesID == namespecies[i])
            j = i
      }
      return(j)
    }
    
    # visibility control for species
    observeEvent(input$speciesID, {
      shinyjs::show("OBJtableSpecies")
      shinyjs::hide("OBJfileRov")
      shinyjs::hide("OBJRov")
      I = setSpecie()
      
      if(natmospecies[I] == 2)
        shinyjs::show("OBJfileRov")
      
      # define print for species
      # output$Inspecies <- renderDataTable({
      #   if(!is.null(input$speciesID)){
      #     aux = setSpecies(I-1)
      #
      #     name = namespecies[I]
      #     natoms = natmospecies[I]
      #     nfreqs = nfreqspecies[I]
      #     symmetry = symmetryspecies[I]
      #     type = typespecies[I]
      #     energy = aux[1,1]
      #     mass = aux[2,1]
      #     x = aux[3,1]
      #     y = aux[4,1]
      #     z = aux[5,1]
      #     frequencies = ""
      #     if(nfreqspecies[I] > 0)
      #       frequencies = aux[6,1]
      #
      #     if(natmospecies[I] > 1){
      #       j = natmospecies[I]
      #       if(natmospecies[I] < nfreqspecies[I])
      #         j = nfreqspecies[I]
      #       for (i in 2:j) {
      #         name = c(name,"")
      #         natoms = c(natoms,"")
      #         nfreqs = c(nfreqs,"")
      #         symmetry = c(symmetry,"")
      #         type = c(type,"")
      #         energy = c(energy,"")
      #         if(natmospecies[I] > nfreqspecies[I]){
      #           mass = c(mass,aux[2,i])
      #           x = c(x,aux[3,i])
      #           y = c(y,aux[4,i])
      #           z = c(z,aux[5,i])
      #           if(i <= nfreqspecies[I])
      #             frequencies = c(frequencies,aux[6,i])
      #           else
      #             frequencies = c(frequencies,"")
      #         }
      #         else{
      #           frequencies = c(frequencies,aux[6,i])
      #           if(i <= natmospecies[I]){
      #             mass = c(mass,aux[2,i])
      #             x = c(x,aux[3,i])
      #             y = c(y,aux[4,i])
      #             z = c(z,aux[5,i])
      #           }
      #           else{
      #             mass = c(mass,"")
      #             x = c(x,"")
      #             y = c(y,"")
      #             z = c(z,"")
      #           }
      #         }
      #       }
      #     }
      #
      #
      #     df = data.frame(name, natoms, nfreqs, symmetry, type, energy, mass, x, y, z, frequencies)
      #     datatable(
      #       df,
      #       options = list(scrollX = TRUE, pageLength = 25),
      #     )
      #   }
      # })
      
      # define print for species
      aux = setSpecies(I-1)
      geo = ""
      for (i in 1:natmospecies[I])
        geo = paste(geo,"<tr><td>",aux[2,i],"&emsp;</td><td>",aux[3,i],"&emsp;</td><td>",aux[4,i],"&emsp;</td><td>",aux[5,i],"</td></tr>")
      freq = ""
      if(nfreqspecies[I] > 0)
        for (i in 1:nfreqspecies[I])
          freq = paste(freq,"<tr><td>",aux[6,i],"</td></tr>")
      
      output$Inspecies2 <- renderText({
        paste0("<h4><b>Name: </b>",namespecies[I],"&emsp;<b>&#8470; Atoms: </b>",natmospecies[I],"&emsp;<b>&#8470; Frequencies: </b>",nfreqspecies[I],
               "&emsp;<b>Symmetry: </b>",symmetryspecies[I],"&emsp;<b>Type: </b>",typespecies[I],"<br><br><b>Energy: </b>",aux[1,1],
               "<br><br><b>Geometry:</b><br>",
               "<table>
                      <tr>
                        <th>mass:</th>
                        <th>x:</th>
                        <th>y:</th>
                        <th>z:</th>
                      </tr>",geo,
               "</table><br>",
               "<table>
                      <tr>
                        <th>Frequencies:</th>
                      </tr>",freq,
               "</table>"
        )
      })
    })
    
    # input for rovibrational levels
    output$RovLevels <- reactive({
      I = setSpecie()
      file <- input$File3
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "txt", "Error: please upload a .txt file for the rovibrational levels"), errorClass = "errorClass")
      
      rovlevels <- read.table(file$datapath, header = FALSE, fill = TRUE, blank.lines.skip = TRUE, strip.white = TRUE, sep = " ", numerals = "no.loss", colClasses = c("character","numeric"))
      
      name <- rovlevels[1,1]
      di <- rovlevels[2:6,2]
      Eref <- rovlevels[7,2]
      De <- rovlevels[8,2]
      Req <- rovlevels[9,2]
      
      if(name == namespecies[I]){
        extractRov(I-1,di,Eref,De,Req)
        
        renames <-function(){
          aux = matrix(0,50,2)
          for (i in 0:49) {
            aux[i+1,1] = paste0("J",i)
            aux[i+1,2] = paste0("v",i)
          }
          return(aux)
        }
        
        # define print for Rvl table
        output$InRvl <- renderDataTable({
          rname = renames()
          datatable(
            format(setRvl(I-1)),
            colnames = rname[1:50,1],
            rownames = rname[1:50,2],
            options = list(scrollX = TRUE, pageLength = 10),
          )
        })
        
        specs = setSpecConst(I-1)
        
        # define print for Spec constants
        output$SpecConst <- renderText({
          paste("<br><b>Y<sub>(0,0)</sub>:</b>",format(specs[1,1]),"&nbsp &nbsp <b>&#969<sub>e</sub>:</b>",format(specs[2,1]),"&nbsp &nbsp <b>&#969<sub>e</sub>X<sub>e</sub>:</b>",format(specs[3,1]),"&nbsp &nbsp <b>&#969<sub>e</sub>Y<sub>e</sub>:</b>",format(specs[4,1]),"&nbsp &nbsp <b>&#969<sub>e</sub>Z<sub>e</sub>:</b>",format(specs[5,1])
                ,"<br><br><b>B<sub>e</sub>:</b>",format(specs[1,2]),"&nbsp &nbsp <b>&#945<sub>e</sub>:</b>",format(specs[2,2]),"&nbsp &nbsp <b>&#947<sub>e</sub>:</b>",format(specs[3,2]),"&nbsp &nbsp <b>&#951<sub>e</sub>:</b>",format(specs[4,2])
                ,"<br><br><b>D<sub>e</sub>:</b>",format(specs[1,3]),"&nbsp &nbsp <b>&#946<sub>e</sub>:</b>",format(specs[2,3]),"&nbsp &nbsp <b>&#916<sub>e</sub>:</b>",format(specs[3,3])
                ,"<br><br><b>F<sub>e</sub>:</b>",format(specs[1,4]),"&nbsp &nbsp <b>Y<sub>(1,3)</sub>:</b>",format(specs[2,4])
                ,"<br><br><b>H<sub>e</sub>:</b>",format(specs[1,5])
                ,"<br><br><b>Y<sub>(0,10)</sub>:</b>",format(specs[1,11]),"&nbsp &nbsp <b>Y<sub>(1,8)</sub>:</b>",format(specs[2,9]),"&nbsp &nbsp <b>Y<sub>(2,6)</sub>:</b>",format(specs[3,7]),"&nbsp &nbsp <b>Y<sub>(3,4)</sub>:</b>",format(specs[4,5]))
        })
        
        shinyjs::show("OBJRov")
        shinyjs::hide("OBJtableSpecies")
      }
      else
        shinyjs::hide("OBJRov")
      
      # define print for specie with Rvl
      aux = setSpecies(I-1)
      geo = ""
      for (i in 1:natmospecies[I])
        geo = paste(geo,"<tr><td>",aux[2,i],"&emsp;</td><td>",aux[3,i],"&emsp;</td><td>",aux[4,i],"&emsp;</td><td>",aux[5,i],"</td></tr>")
      freq = ""
      if(nfreqspecies[I] > 0)
        for (i in 1:nfreqspecies[I])
          freq = paste(freq,"<tr><td>",aux[6,i],"</td></tr>")
      
      output$Inspecie <- renderText({
        paste0("<h4><b>Name: </b>",namespecies[I],"&emsp;<b>&#8470; Atoms: </b>",natmospecies[I],"&emsp;<b>&#8470; Frequencies: </b>",nfreqspecies[I],
               "&emsp;<b>Symmetry: </b>",symmetryspecies[I],"&emsp;<b>Type: </b>",typespecies[I],"<br><br><b>Energy: </b>",aux[1,1],
               "<br><br><b>Geometry:</b><br>",
               "<table>
                      <tr>
                        <th>mass:</th>
                        <th>x:</th>
                        <th>y:</th>
                        <th>z:</th>
                      </tr>",geo,
               "</table><br>",
               "<table>
                      <tr>
                        <th>Frequencies:</th>
                      </tr>",freq,
               "</table>"
        )
      })
      
      # output$Inspecie <- renderDataTable({
      #   if(!is.null(input$speciesID)){
      #     i = 1
      #     for (k in 1:length(species[,1])) {
      #       if(!is.na(species[k,1]))
      #         if(namespecies[I] == species[k,1])
      #           i=k
      #     }
      #     if(species[i,2] > species[i,3])
      #       j = i + species[i,2] - 1
      #     else
      #       j = i + species[i,3] - 1
      #
      #     # define print for new energy
      #     output$NewEnergy <- renderText({
      #       paste("<h4>Current energy:",species[i,6])
      #     })
      #
      #     datatable(
      #       species[i:j,1:11],
      #       rownames = FALSE,
      #       options = list(scrollX = TRUE, pageLength = 25),
      #     )
      #   }
      # })
      
      validate(need(ext != "txt", "The .txt file for rovibrational levels has been read"))
    })
    
    # add the (v,J) level
    observeEvent(input$EnterRov, {
      I = setSpecie()
      aux = NULL
      
      energya = addRvl(I-1, 0, 0)
      energyb = addRvl(I-1, input$levelv, input$levelJ)
      aux1 = paste("<h4>Original energy:",energya,"<br><br>New energy:<font color=red>",energyb)
      
      # define print for new energy
      output$NewEnergy <- renderText({
        aux1
      })
      
      # define print for specie with Rvl
      aux = setSpecies(I-1)
      geo = ""
      for (i in 1:natmospecies[I])
        geo = paste(geo,"<tr><td>",aux[2,i],"&emsp;</td><td>",aux[3,i],"&emsp;</td><td>",aux[4,i],"&emsp;</td><td>",aux[5,i],"</td></tr>")
      freq = ""
      if(nfreqspecies[I] > 0)
        for (i in 1:nfreqspecies[I])
          freq = paste(freq,"<tr><td>",aux[6,i],"</td></tr>")
      
      output$Inspecie <- renderText({
        paste0("<h4><b>Name: </b>",namespecies[I],"&emsp;<b>&#8470; Atoms: </b>",natmospecies[I],"&emsp;<b>&#8470; Frequencies: </b>",nfreqspecies[I],
               "&emsp;<b>Symmetry: </b>",symmetryspecies[I],"&emsp;<b>Type: </b>",typespecies[I],"<br><br><b>Energy: </b>",aux[1,1],
               "<br><br><b>Geometry:</b><br>",
               "<table>
                      <tr>
                        <th>mass:</th>
                        <th>x:</th>
                        <th>y:</th>
                        <th>z:</th>
                      </tr>",geo,
               "</table><br>",
               "<table>
                      <tr>
                        <th>Frequencies:</th>
                      </tr>",freq,
               "</table>"
        )
      })
    })
    
    therespecies <<- TRUE
    validate(need(ext != "ods", "The .ods file for species has been read"))
  })
  
  # input for reactions
  output$Reactions <- reactive({
    therereactions <<- FALSE
    shinyjs::hide("OBJoutput")
    shinyjs::hide("OBJLets")
    file <- input$File2
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "ods", "Error: please upload a .ods file for the reactions"), errorClass = "errorClass")
    
    reactions <<- read_ods(path =  file$datapath, col_names = TRUE)
    
    # define print for reactions
    output$Inreactions <- renderDataTable({
      datatable(
        reactions,
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 25),
      )
    })
    
    nreactions <<- length(reactions$reactant1[!is.na(reactions$reactant1)])
    react1 <<- reactions$reactant1[!is.na(reactions$reactant1)]
    react2 <<- reactions$reactant2[!is.na(reactions$reactant2)]
    tss <<- reactions$ts[!is.na(reactions$ts)]
    prod1 <<- reactions$product1[!is.na(reactions$product1)]
    prod2 <<- reactions$product2[!is.na(reactions$product2)]
    psi <- reactions$psi[!is.na(reactions$psi)]
    
    extractReactions(nreactions,react1,react2,tss,prod1,prod2,psi)
    
    therereactions <<- TRUE
    validate(need(ext != "ods", "The .ods file for reactions has been read"))
  })
  
  # visibility control for file1
  observeEvent(input$File1, {
    shinyjs::show("OBJspecies")
    if(therereactions == TRUE)
      shinyjs::show("OBJreactions")
    shinyjs::hide("OBJinputs")
    shinyjs::hide("OBJerrorinputs")
    shinyjs::hide("OBJcalculate")
    shinyjs::hide("OBJerrorcalc")
  })
  
  # visibility control for file2
  observeEvent(input$File2, {
    if (therespecies == TRUE)
      shinyjs::show("OBJspecies")
    shinyjs::show("OBJreactions")
    shinyjs::hide("OBJinputs")
    shinyjs::hide("OBJerrorinputs")
    shinyjs::hide("OBJcalculate")
    shinyjs::hide("OBJerrorcalc")
  })
  
  # visibility control for showinputs
  observeEvent(input$Showinputs, {
    shinyjs::hide("OBJspecies")
    shinyjs::hide("OBJreactions")
    shinyjs::hide("OBJcalculate")
    shinyjs::hide("OBJerrorcalc")
    shinyjs::hide("OBJLets")
    if (therespecies == TRUE & therereactions == TRUE)
      shinyjs::show("OBJinputs")
    else
      shinyjs::show("OBJerrorinputs")
  })
  
  # error message for calculation
  output$ErrorCalc <- renderText({
    "Error: species names in input files must be the same"
  })
  
  # visibility control for calculate
  observeEvent(input$Calculate, {
    shinyjs::hide("OBJspecies")
    shinyjs::hide("OBJreactions")
    shinyjs::hide("OBJinputs")
    shinyjs::hide("OBJerrorinputs")
    shinyjs::hide("OBJLets")
    
    if (therespecies == TRUE & therereactions == TRUE){
      
      # verify if there's some errors in species names
      count = 0
      for(i in 1:nreactions){
        if(react2[i] == "ZERO")
          count = count + 1
        if(prod2[i] == "ZERO")
          count = count + 1
        for (j in 1:nspecies) {
          if(react1[i] == namespecies[j])
            count = count + 1
          if(react2[i] == namespecies[j])
            count = count + 1
          if(tss[i] == namespecies[j])
            count = count + 1
          if(prod1[i] == namespecies[j])
            count = count + 1
          if(prod2[i] == namespecies[j])
            count = count + 1
        }
      }
      
      if(count == 5*nreactions){
        calculateRate()
        
        # define which reaction to show
        output$ReactionID <- renderUI({
          selectInput("reactionID", "Choose reaction by TS:", tss, width = '100%')
        })
        
        # set the actual TS to show
        setTSi <- function(){
          j = 1
          if(is.null(input$reactionID) == FALSE)
            for (i in 1:nreactions)
              if(input$reactionID == tss[i])
                j = i
          return(j)
        }
        
        # define the selected reaction
        output$ReactionPrint <- renderUI({
          I = setTSi()
          tags$h4(strong("Reaction"),strong(I),strong(":"),react1[I],"+",react2[I],icon("long-arrow-alt-right"),tss[I],icon("long-arrow-alt-right"),prod1[I],"+",prod2[I])
        })
        
        # set the reaction rate graph
        plotRateInput = function(){
          I = setTSi()
          auxrate = setRate(I-1)
          par(mar=c(5,6,4,1)+.1)
          if (!input$checkNa){
            laby = expression(log(k(cm^{3}~mol^{-1}~s^{-1})))
            if(react2[I] == "ZERO")
              laby = expression(log(k(s^{-1})))
            matplot(auxrate[,1],auxrate[,2:5], type = "l", col = 2:5, xlab = expression(10000/T(K^{-1})), ylab = laby)
            legend("topright",
                   inset = 0.05,
                   legend = c("k(T)", expression(k[W](T)), expression(k[E](T)), expression(k[SCT](T))),
                   lty = 1:4,
                   col = c(2, 3, 4, 5),
                   lwd = 2)
          }
          else{
            auxrate[,2:5] = 10^(auxrate[,2:5])
            auxrate[,2:5] = auxrate[,2:5]/6.022e23
            auxrate[,2:5] = log10(auxrate[,2:5])
            laby = expression(log(k(cm^{3}~s^{-1})))
            if(react2[I] == "ZERO")
              laby = expression(log(k(mol~s^{-1})))
            matplot(auxrate[,1],auxrate[,2:5], type = "l", col = 2:5, xlab = expression(10000/T(K^{-1})), ylab = laby)
            legend("topright",
                   inset = 0.05,
                   legend = c("k(T)", expression(k[W](T)), expression(k[E](T)), expression(k[SCT](T))),
                   lty = 1:4,
                   col = c(2, 3, 4, 5),
                   lwd = 2)
          }
        }
        
        output$PlotRate <- renderPlot({plotRateInput()
        })
        
        # define print for arrhenius
        output$Arrh <- renderUI({
          I = setTSi()
          arrh = setArrhenius(I-1)
          if(react2[I] == "ZERO"){
            if (!input$checkNa){
              withMathJax(
                sprintf('$$Arrhenius \\ form: \\ (s^{-1})$$'),
                sprintf("$$k(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",arrh[1,1]/10^floor(log10(arrh[1,1])),floor(log10(arrh[1,1])),arrh[1,2],-arrh[1,3]/1000),
                sprintf("$$k_W(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",arrh[2,1]/10^floor(log10(arrh[2,1])),floor(log10(arrh[2,1])),arrh[2,2],-arrh[2,3]/1000),
                sprintf("$$k_E(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",arrh[3,1]/10^floor(log10(arrh[3,1])),floor(log10(arrh[3,1])),arrh[3,2],-arrh[3,3]/1000),
                sprintf("$$k_{SCT}(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",arrh[4,1]/10^floor(log10(arrh[4,1])),floor(log10(arrh[4,1])),arrh[4,2],-arrh[4,3]/1000)
              )
            }
            else{
              withMathJax(
                sprintf('$$Arrhenius \\ form: \\ (mol \\ s^{-1})$$'),
                sprintf("$$k(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",(arrh[1,1]/6.022e23)/10^floor(log10((arrh[1,1]/6.022e23))),floor(log10((arrh[1,1]/6.022e23))),arrh[1,2],-arrh[1,3]/1000),
                sprintf("$$k_W(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",(arrh[2,1]/6.022e23)/10^floor(log10((arrh[2,1]/6.022e23))),floor(log10((arrh[2,1]/6.022e23))),arrh[2,2],-arrh[2,3]/1000),
                sprintf("$$k_E(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",(arrh[3,1]/6.022e23)/10^floor(log10((arrh[3,1]/6.022e23))),floor(log10((arrh[3,1]/6.022e23))),arrh[3,2],-arrh[3,3]/1000),
                sprintf("$$k_{SCT}(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",(arrh[4,1]/6.022e23)/10^floor(log10((arrh[4,1]/6.022e23))),floor(log10((arrh[4,1]/6.022e23))),arrh[4,2],-arrh[4,3]/1000)
              )
            }
          }
          else{
            if (!input$checkNa){
              withMathJax(
                sprintf('$$Arrhenius \\ form: \\ (cm^3 \\ mol^{-1} \\ s^{-1})$$'),
                sprintf("$$k(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",arrh[1,1]/10^floor(log10(arrh[1,1])),floor(log10(arrh[1,1])),arrh[1,2],-arrh[1,3]/1000),
                sprintf("$$k_W(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",arrh[2,1]/10^floor(log10(arrh[2,1])),floor(log10(arrh[2,1])),arrh[2,2],-arrh[2,3]/1000),
                sprintf("$$k_E(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",arrh[3,1]/10^floor(log10(arrh[3,1])),floor(log10(arrh[3,1])),arrh[3,2],-arrh[3,3]/1000),
                sprintf("$$k_{SCT}(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",arrh[4,1]/10^floor(log10(arrh[4,1])),floor(log10(arrh[4,1])),arrh[4,2],-arrh[4,3]/1000)
              )
            }
            else{
              withMathJax(
                sprintf('$$Arrhenius \\ form: \\ (cm^3 \\ s^{-1})$$'),
                sprintf("$$k(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",(arrh[1,1]/6.022e23)/10^floor(log10((arrh[1,1]/6.022e23))),floor(log10((arrh[1,1]/6.022e23))),arrh[1,2],-arrh[1,3]/1000),
                sprintf("$$k_W(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",(arrh[2,1]/6.022e23)/10^floor(log10((arrh[2,1]/6.022e23))),floor(log10((arrh[2,1]/6.022e23))),arrh[2,2],-arrh[2,3]/1000),
                sprintf("$$k_E(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",(arrh[3,1]/6.022e23)/10^floor(log10((arrh[3,1]/6.022e23))),floor(log10((arrh[3,1]/6.022e23))),arrh[3,2],-arrh[3,3]/1000),
                sprintf("$$k_{SCT}(T) = %.2f \\times 10^{%d} \\ T^{%.3f}exp(%.1f \\ kcal \\ mol^{-1}/RT)$$",(arrh[4,1]/6.022e23)/10^floor(log10((arrh[4,1]/6.022e23))),floor(log10((arrh[4,1]/6.022e23))),arrh[4,2],-arrh[4,3]/1000)
              )
            }
          }
        })
        
        # set the MEP graph
        plotMEPInput = function(){
          I = setTSi()
          aux = setMEP(I-1)
          matplot(aux[,1],aux[,2:3], type = "l", col = 2:3, xlab = "reaction coordinate (a.u.)", ylab = expression(kcal~mol^{-1}))
          legend("topright",
                 inset = 0.05,
                 legend = c(expression(V[MEP](s)), expression(V[a]^{G}*(s))),
                 lty = 1:2,
                 col = c(2, 3),
                 lwd = 2)
        }
        
        output$PlotMEP <- renderPlot({plotMEPInput()
        })
        
        #set which TP to show
        setWhichTP = function(j){
          I = setTSi()
          reac1 = reac2 = pro1 = pro2 = matrix(0, 60, 4)
          for (i in 1:nspecies) {
            if(namespecies[i] == react1[I])
              reac1 = setThermProp(i-1)
            if(namespecies[i] == react2[I])
              reac2 = setThermProp(i-1)
            if(namespecies[i] == prod1[I])
              pro1 = setThermProp(i-1)
            if(namespecies[i] == prod2[I])
              pro2 = setThermProp(i-1)
          }
          aux = matrix(0, 60, 5)
          aux[,1] = reac1[,1]
          aux[,2] = reac1[,j]
          aux[,3] = reac2[,j]
          aux[,4] = pro1[,j]
          aux[,5] = pro2[,j]
          return(aux)
        }
        
        # set the Enthalpy (Eao) equation
        output$EqEao <- renderUI({
          withMathJax(
            sprintf('$$E=k_BT \\left(\\frac{\\partial \\ \\ln \\ Q}{\\partial \\ \\ln \\ T}\\right)_\\upsilon$$'),
            sprintf('$$H = E+pV $$'),
            sprintf('$$H = R \\left(a_6 + \\sum_{i=1}^{5}\\frac{a_i T^i}{i}\\right)$$'),
            sprintf('$$H = R \\left(a_8 + a_2 \\ln T + \\sum_{i=1,~i \\neq 2}^{7}\\frac{a_i T^{i-2}}{i-2}\\right)$$')
          )
        })
        
        # set the Enthalpy (Eao) graph
        plotEaoInput = function(){
          I = setTSi()
          aux = setWhichTP(2)
          matplot(aux[,1],aux[,2:5], type = "l", col = 2:5, xlab = "T(K)", ylab = expression(cal~mol^{-1}))
          legend("topright",
                 inset = 0.05,
                 legend = c(react1[I],react2[I],prod1[I],prod2[I]),
                 lty = 1:4,
                 col = c(2, 3, 4, 5),
                 lwd = 2)
        }
        
        output$PlotEao <- renderPlot({plotEaoInput()
        })
        
        # set the Entropy (Sto) equation
        output$EqSto <- renderUI({
          withMathJax(
            sprintf('$$S=k_B \\ \\ln \\ Q + k_B \\left(\\frac{\\partial \\ \\ln \\ Q}{\\partial \\ \\ln \\ T}\\right)_\\upsilon$$'),
            sprintf('$$S = R \\left(a_1 \\ln T + a_7 + \\sum_{i=2}^{5}\\frac{a_i T^{i-1}}{i-1}\\right)$$'),
            sprintf('$$S = R \\left(a_3 \\ln T + a_9 + \\sum_{i=1,~i \\neq 3}^{7}\\frac{a_i T^{i-3}}{i-3}\\right)$$')
          )
        })
        
        # set the Entropy (Sto) graph
        plotStoInput = function(){
          I = setTSi()
          aux = setWhichTP(3)
          matplot(aux[,1],aux[,2:5], type = "l", col = 2:5, xlab = "T(K)", ylab = expression(cal~K^{-1}~mol^{-1}))
          legend("topright",
                 inset = 0.05,
                 legend = c(react1[I],react2[I],prod1[I],prod2[I]),
                 lty = 1:4,
                 col = c(2, 3, 4, 5),
                 lwd = 2)
        }
        
        output$PlotSto <- renderPlot({plotStoInput()
        })
        
        # set the Heat Capacity (Cpo) equation
        output$EqCpo <- renderUI({
          withMathJax(
            sprintf('$$c_\\upsilon=k_B \\left(\\frac{\\partial \\ \\ln \\ Q}{\\partial \\ \\ln \\ T}\\right)_\\upsilon + k_B \\left(\\frac{\\partial^2 \\ \\ln \\ Q}{\\partial (\\ln \\ T)^2}\\right)_\\upsilon$$'),
            sprintf('$$c_p = c_\\upsilon + R $$'),
            sprintf('$$c_p = R \\sum_{i=1}^{5} a_i T^{i-1}$$'),
            sprintf('$$c_p = R \\sum_{i=1}^{7} a_i T^{i-3}$$')
          )
        })
        
        # set the Heat Capacity (Cpo) graph
        plotCpoInput = function(){
          I = setTSi()
          aux = setWhichTP(4)
          matplot(aux[,1],aux[,2:5], type = "l", col = 2:5, xlab = "T(K)", ylab = expression(cal~K^{-1}~mol^{-1}))
          legend("topright",
                 inset = 0.05,
                 legend = c(react1[I],react2[I],prod1[I],prod2[I]),
                 lty = 1:4,
                 col = c(2, 3, 4, 5),
                 lwd = 2)
        }
        
        output$PlotCpo <- renderPlot({plotCpoInput()
        })
        
        # set download button
        output$DownloadData <- downloadHandler(
          filename = function(){
            paste0("alldata.zip")},
          content = function(file){
            #go to a temp dir to avoid permission issues
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            files <- NULL;
            
            #loop through the sheets
            for (i in 1:nreactions){
              #write each sheet of reaction rate to a .ods file, save the name
              fileName <- paste("rate",tss[i],".ods",sep = "")
              rate <- setRate(i-1)
              df <- data.frame(`T` = rate[,1], `k(T)` = rate[,2], `kW(T)` = rate[,3], `kE(T)` = rate[,4], `kSCT(T)` = rate[,5], check.names = FALSE)
              write_ods(df, fileName, row_names = FALSE, col_names = TRUE)
              files <- c(fileName,files)
              
              #write each sheet of MEP to a .ods file, save the name
              fileName1 <- paste("mep",tss[i],".ods",sep = "")
              mep <- setMEP(i-1)
              df <- data.frame(s = mep[,1], `Vmep(s)` = mep[,2], `VaG(s)` = mep[,3], check.names = FALSE)
              write_ods(df, fileName1, row_names = FALSE, col_names = TRUE)
              files <- c(fileName1,files)
            }
            for (i in 1:nspecies) {
              #write each sheet of thermodynamic properties to a .ods file, save the name
              fileName <- paste("thermprop",namespecies[i],".ods",sep = "")
              tprop <- setThermProp(i-1)
              df <- data.frame(`T` = tprop[,1], `Enthalpy(T)` = tprop[,2], `Entropy(T)` = tprop[,3], `Heat capacity(T)` = tprop[,4], check.names = FALSE)
              write_ods(df, fileName, sheet = "Thermodynamic Properties", row_names = FALSE, col_names = TRUE)
              df <- data.frame(tprop[1:7,6:10])
              write_ods(df, fileName, sheet = "NASA Coefficients", append = TRUE, row_names = FALSE, col_names = FALSE)
              files <- c(fileName,files)
            }
            #create the zip file
            zip(file,files)
          }
        )
        
        shinyjs::show("OBJcalculate")
        shinyjs::show("OBJoutput")
      }
      else
        shinyjs::show("OBJerrorcalc")
    }
    else{
      shinyjs::show("OBJerrorinputs")
      shinyjs::hide("OBJoutput")
    }
  })
}