#for distribution use
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
if(!require(shinyFiles)){
  install.packages("shinyFiles")
  library(shinyFiles)
}
if(!require(vroom)){
  install.packages("vroom")
  library(vroom)
}
if(!require(shinycssloaders)){
  install.packages("shinycssloaders")
  library(shinycssloaders)
}

#for basic logic use
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(UpSetR)){
  install.packages("UpSetR")
  library(UpSetR)
}
if(!require(VennDiagram)){
  install.packages("VennDiagram")
  library(VennDiagram)
}


###################################################################################
#================================     UI      ====================================#
###################################################################################

ui <- fluidPage(
  
 
  #navbar icon pic
  navbarPage(title = "PPI Tool",
    #title=div(img(src="MyIcon.png",height = 25, width = 30), "PPI Tool"),
          
            
             
             #background pic
             tags$style(
               'div[data-value="Home"]{ 
                 width: 100%; height: 950px; 
                 background-image: url("background_img_mainpage.png");
               }'
             ),
             
             #============================== Main page ================================#  
             tabPanel("Home",
                      
                      
                      titlePanel(h1("PPI tool",align="center",style="color: skyblue")), 
                      
                      titlePanel(h4("This tool is for protein-protein interactor screening from mass spectrometry data. ",align="center",style="color: orange")),
                      br(),
                      titlePanel(h4("It can be applied to both affinity-purification mass spectrometry data and",align="center")),    
                      titlePanel(h4("proximity-labeling mass spectrometry data. PSM information is required but only ",align="center")), 
                      titlePanel(h4("label-free DDA data has been tested.",align="center")), 
                      br(),
                      titlePanel(h4("PPI tool started by modifying Fisher’s exact test. We suppose the distribution ",align="center")),
                      titlePanel(h4("of true and false interactors from PPI studies follows the same ",align="center")),
                      titlePanel(h4("hypergeometric distribution. Generally, data with biological replicates are ",align="center")),
                      titlePanel(h4("recommended but data without replicates are also suitable.",align="center")),
                   
                       
             ), 
             
             
             #============================== Analysis page ================================#    
             tabPanel("Analysis",
                      
                      fluidRow(
                        
                        sidebarLayout(
                          sidebarPanel(
                            
                            
                            #input data
                            fileInput("file1", "Choose CSV file", accept = ".csv"),
                            
                            # Number of repeats
                            sliderInput("repeats", "Number of replicates", value = 3, min = 1, max = 5),
                            
                            #Fisher p-value
                            numericInput("pvalue", "1. P-value for Fisher cut-off", value= 0.05),
                            
                            #number of comparing set
                            numericInput("num1", "2. The bait vs. # of background group for Fisher.test (Set 0 for the overview of all significant lists in compare to all backgrounds and Upset plot)", value = 0, min = 0, max = 999),
                            
                            actionButton("goButton", "Run"),
                            
                            
                          ),
                          
                          
                            
                            mainPanel(
                              
                              tabsetPanel(
                                
                                #loading effect withspinner()
                                tabPanel("List of prey",withSpinner(tableOutput("contents")),downloadButton('tableDown',label = 'Download Table')),
                                
                                tabPanel("Upset plot",withSpinner(plotOutput("contents2",width = "100%",height = "800px")),downloadButton('upsetplotDown',label = 'Download Plot')),
                                
                                tabPanel("Interactions chart",withSpinner(tableOutput("contents3")),downloadButton('tableDown_intersect',label = 'Download Table')),
                                
                              )
                            ) 
                            
                          
                          
                        )
                        
                      )
             ),
             
             #============================== Documentation ================================#
             
             tabPanel("Documentation",
                      mainPanel(
                        # img(src="FAFUicon.jpg", height = 100, width = 100),
                        br(),
                        
                        
                        h2("Preparation of data sheet"),
                        br(),
                        h4("Name your .csv file as your bait protein name, e.g. 'APC1.csv'. "),
                        
                        
                        br(),
                        p("Column 1-5: Information about the proteins. At column 2, put single protein ID, the rest column can be any information."),
                        br(),
                        p("Starting from column 6, list the spectral count of bait & background protein."),
                        br(),
                        p("Following the columns, list the data of other baits or backgrounds (same replicates number). "),
                        br(),
                        img(src="doc_pic1.jpg", height = 650, width = 1400),
                        
                        br(),
                        br(),
                        
                        h2("Input Parameters "),
                        br(),
                        
                        p("1. P-value for Fisher cut-off"),
                        br(),
                        p("2. The bait vs. # of background group for '1 on 1' Fisher.test " ),
                        p("   and shows p-value of each proteins  (The following instance uses APC1.csv in input_example folder )" ),
                        br(),
                        img(src="doc_pic2.png"),
                        
                        
                        
                      )
             ),
             
             #============================== About ===============================#
             
             tabPanel("About",
                      mainPanel(
                        
                        
                        
                        p(""),
                        br(),
                        p("The PPI tool was developed by Qiyao Wu under the guidance of Dr. Xueyang Zhang at the proteomics center of Haixia Institute of Science and Technology at Fujian Agriculture and Forestry University."),
                        br(),
                        p("Thanks for the funding support provided by dear Prof. Dr. Chentao Lin. "),
                        br(),
                        p("Please contact Qiyao Wu (260702397@qq.com) for technical issues. "),
                        br(),
                        p("Or Dr. Xueyang Zhang (xueyangzhang2018@163.com) for other stuff. "),
                        br(),
                        br(),
                        img(src="FAFUicon.jpg", height = 100, width = 100),
                        
                        
                        
                        
                      )
             )
             
  )
)

#################################################################################################
#=======================================   Server   ============================================#
#################################################################################################

server <- function(input, output) {
  
  getData1 <- reactive({
    
    input$goButton
    
  
     
    if(input$goButton!= 0){
     
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "csv", "Please upload a csv file"))
      
      #------------------Data processing start------------------#
      
      data1 <- read.csv(file$datapath, header = T)
      
      bait_target_name <- gsub(".csv","",file)
      
      col_number <- ncol(data1)
      
      #take column of data, starting from column 6, 1-5 columns are char. info
      data_num <- data.frame(data1[,c(6:col_number)])
      
      #take bait data
      bait <- data_num[,c(1:input$repeats)]
      
      #take background data
      back <- data_num[,c((input$repeats +1 ):(col_number-5))]
      
      #count the number of background groups
      n <- ncol(back)/input$repeats
      
      #create an empty list for collecting significance of different groups of comparison
      result_cutoff_protein_list <- list()
      
      
      #===================================================================================================#
      #when num1 is not 0, conduct 1 bait vs 1 background job
      
      if(input$num1 != 0) {
        
        i <- input$num1  
        
        
        #According to 'input$repeats', take the specific background group '#i' 
        if (input$repeats == 1 ) {
          back_inloop <- back[,c( (i-1)*1 +1 
                                  
          )
          ]
        } else
          if (input$repeats == 2 ) {
            back_inloop <- back[,c( (i-1)*2 +1, 
                                    (i-1)*2 +2 
                                    
            )
            ]
          } else
            if(input$repeats == 3 ) {
              
              back_inloop <- back[,c( (i-1)*3 +1, 
                                      (i-1)*3 +2, 
                                      (i-1)*3 +3
              )
              ]
            } else 
              if (input$repeats == 4 ) {
                back_inloop <- back[,c( (i-1)*4 +1, 
                                        (i-1)*4 +2, 
                                        (i-1)*4 +3,
                                        (i-1)*4 +4  
                )
                ]
              } else
                if (input$repeats == 5 ) {
                  back_inloop <- back[,c( (i-1)*5 +1, 
                                          (i-1)*5 +2, 
                                          (i-1)*5 +3,
                                          (i-1)*5 +4,
                                          (i-1)*5 +5,
                  )
                  ]
                }
        
        
        df <-cbind(data1[,c(1:5)],
                   bait, 
                   back_inloop
        )
        
        
        #--------------bait background sheet for Fisher_Test--------------#
        
        #take char. info
        in_df <- df[,c(1:5)]
        
        #Total count of background groups
        background_total <- sum(back_inloop)
        
        #Total count of bait groups
        bait_total <- sum(bait)
        
        
        #comparing total of backg. and bait, deciding the next process
        if (background_total < bait_total)
          
        {
          back_inloop2 <- back_inloop
          bait_inloop2 <- bait
          
          
          
          background_total_loop2 <- background_total
          bait_total_loop2       <- bait_total
          
        } 
        
        #If bait < backg, start imputing
        if (background_total > bait_total)
          
        {
          
          df2 <- cbind(in_df,bait)
          df2_original <- df2      #For comparing before and after imputation
          
          row_num <- nrow(df2)
          
          rows_of_zeros <- 0
          
          #count all rows of '0s'
          for (x in 1:row_num) 
            
          {
            matrix1 <- as.matrix(df2[x,c(6:(6+input$repeats-1))])
            
            
            
            if ( sum(matrix1) ==0 ) {
              rows_of_zeros = rows_of_zeros +1
            }
            
          }
          
          #fold of bait and backg. 
          fold1 <- (background_total - bait_total)/rows_of_zeros
          
          #define imputing number (at maximum of 4)
          if (fold1 < 1)
          {impute_num <- 1
          impute_star <- "*"} else 
            if (fold1 > 1 & fold1 < 2 )
            {impute_num <- 2
            impute_star <- "**"} else 
              if (fold1 > 2 & fold1 < 3 )
              {impute_num <- 3
              impute_star <- "***"} else 
                if (fold1 > 3 & fold1 < 4 )
                {impute_num <- 4
                impute_star <- "****"} else 
                  if (fold1 > 4 )
                  {impute_num <- 4
                  impute_star <- "*****"
                  }
          
          
          for (x in 1:row_num) 
          {
            matrix1 <- as.matrix(df2[x,c(6:(6+input$repeats-1))])
            
            if ( sum(matrix1) ==0 )
              
            { 
              matrix1[1,1] <- impute_num      
              df2[x,c(6:8)]    <-    matrix1  
              df2[x,2] <- paste(df2[x,2],sep = "")   
            }
          }
          
            
          in_df <- df2[,c(1:5)] 
          
          back_inloop2 <- back_inloop
          bait_inloop2 <- df2[,c(6:(6+input$repeats-1))]
          
          
          
          background_total_loop2 <- background_total
          bait_total_loop2       <- sum(df2[,c(6:(6+input$repeats-1))]) 
          
        }
        
        
        
        
        
        
        #row background sum
        in_df$background <- apply(back_inloop2,1,sum)
        
        
        #background othesr
        in_df$background_others <- background_total_loop2 - in_df$background
        
        #bait row sum  
        in_df$baited <-  apply(bait_inloop2,1,sum)
        
        #bait 
        in_df$bait <- in_df$baited -  in_df$background
        
        
        in_df <- subset(in_df,select = -c(baited))
        
        #negative to zero
        in_df$bait[in_df$bait < 0 ] = 0
        
        
        #bait_others 
        in_df$bait_others <- bait_total_loop2 - background_total_loop2 - in_df$bait
        
        #negative to zero
        in_df$bait_others[in_df$bait_others < 0 ] = 0
        
        
        
        #rename col
        colnames(in_df)[6:9] <- c(
          paste(gsub('.{2}$','',colnames(back_inloop2)[1]),"background",sep = "_"),
          paste(gsub('.{2}$','',colnames(back_inloop2)[1]),"background_others",sep = "_"),
          paste(gsub('.{2}$','',colnames(bait_inloop2)[1]),"bait",sep = "_"),
          paste(gsub('.{2}$','',colnames(bait_inloop2)[1]),"bait_others",sep = "_")      
        )
        
        
        #===================================#
        #---------start Fisher loop---------#
        #===================================#
        
        
        data3 <- in_df[,6:9]
        
        
        k_row <- nrow(data3)
        
        
        Fisher_P <-data.frame()
        
        
        for (k in 1:k_row) {
          
          #~~~#
          tryCatch({  
            #~~~#
            
            #构建fisher矩阵
            datainput <- matrix(as.numeric(c(data3[k,1],data3[k,2],data3[k,3],data3[k,4])), 2,
                                dimnames = list(c("protein","others"),
                                                c("background","bait")
                                )
            )
            #Fisher Test
            dataoutput <- fisher.test(datainput, alternative = "less")
            
            #collect P
            Fisher_P[k,1] <- dataoutput$p.value
            
            #~~~#  #
          }, error=function(e){print(paste("=======","Error occurred in",result_out[i],data2[k,1],"=======",sep = " "))})  
          #~~~# 
        }
        
        colnames(Fisher_P) <- "Fisher_P"
        
        result_cutoff_protein_list <- data.frame(cbind(in_df,Fisher_P))
        
        
        
        result_cutoff_protein_list
        
        
      } else
        if (input$num1 == 0) {
          
          #===================================================================================================================#
          #When num1 = 0 
          for (i in 1:n) {
            
            
            
            if (input$repeats == 1 ) {
              back_inloop <- back[,c( (i-1)*1 +1 
                                      
              )
              ]
            } else
              if (input$repeats == 2 ) {
                back_inloop <- back[,c( (i-1)*2 +1, 
                                        (i-1)*2 +2 
                                        
                )
                ]
              } else
                if(input$repeats == 3 ) {
                  
                  back_inloop <- back[,c( (i-1)*3 +1, 
                                          (i-1)*3 +2, 
                                          (i-1)*3 +3
                  )
                  ]
                } else 
                  if (input$repeats == 4 ) {
                    back_inloop <- back[,c( (i-1)*4 +1, 
                                            (i-1)*4 +2, 
                                            (i-1)*4 +3,
                                            (i-1)*4 +4  
                    )
                    ]
                  } else
                    if (input$repeats == 5 ) {
                      back_inloop <- back[,c( (i-1)*5 +1, 
                                              (i-1)*5 +2, 
                                              (i-1)*5 +3,
                                              (i-1)*5 +4,
                                              (i-1)*5 +5,
                      )
                      ]
                    }
            
            
            df <-cbind(data1[,c(1:5)],
                       bait, 
                       back_inloop
            )
            
            
            #--------------bait background sheet for Fisher_Test--------------#
            
            
            in_df <- df[,c(1:5)]
            
            
            background_total <- sum(back_inloop)
            
            
            bait_total <- sum(bait)
            
            
            
            if (background_total < bait_total)
              
            {
              back_inloop2 <- back_inloop
              bait_inloop2 <- bait
              
              
              
              background_total_loop2 <- background_total
              bait_total_loop2       <- bait_total
              
            } 
            
            
            if (background_total > bait_total)
              
            {
              
              df2 <- cbind(in_df,bait)
              df2_original <- df2     
              
              row_num <- nrow(df2)
              
              rows_of_zeros <- 0
              
              
              for (x in 1:row_num) 
                
              {
                matrix1 <- as.matrix(df2[x,c(6:(6+input$repeats-1))])
                
                
                
                if (sum (matrix1) ==0 ) {
                  rows_of_zeros = rows_of_zeros +1
                }
                
              }
              
              
              fold1 <- (background_total - bait_total)/rows_of_zeros
              
              
              if (fold1 < 1)
              {impute_num <- 1
              impute_star <- "*"} else 
                if (fold1 > 1 & fold1 < 2 )
                {impute_num <- 2
                impute_star <- "**"} else 
                  if (fold1 > 2 & fold1 < 3 )
                  {impute_num <- 3
                  impute_star <- "***"} else 
                    if (fold1 > 3 & fold1 < 4 )
                    {impute_num <- 4
                    impute_star <- "****"} else 
                      if (fold1 > 4 )
                      {impute_num <- 4
                      impute_star <- "*****"
                      }
              
              
              for (x in 1:row_num) 
              {
                matrix1 <- as.matrix(df2[x,c(6:(6+input$repeats-1))])
                
                if (sum(matrix1)==0 )
                  
                { 
                  matrix1[1,1] <- impute_num      
                  df2[x,c(6:(6+input$repeats-1))]    <-    matrix1  
                  df2[x,2] <- paste(df2[x,2],impute_star,sep = "")  
                }
              }
              
                
              in_df <- df2[,c(1:5)] 
              
              back_inloop2 <- back_inloop
              bait_inloop2 <- df2[,c(6:(6+input$repeats-1))]
              
              
              
              background_total_loop2 <- background_total
              bait_total_loop2       <- sum(df2[,c(6:(6+input$repeats-1))]) 
              
            }
            
            
            
            
            
            
            
            in_df$background <- apply(back_inloop2,1,sum)
            
            
            
            in_df$background_others <- background_total_loop2 - in_df$background
            
             
            in_df$baited <-  apply(bait_inloop2,1,sum)
            
            
            in_df$bait <- in_df$baited -  in_df$background
            
            
            in_df <- subset(in_df,select = -c(baited))
            
            
            in_df$bait[in_df$bait < 0 ] = 0
            
            
            
            in_df$bait_others <- bait_total_loop2 - background_total_loop2 - in_df$bait
            
            
            in_df$bait_others[in_df$bait_others < 0 ] = 0
            
            
            
            
            colnames(in_df)[6:9] <- c(
              paste(gsub('.{2}$','',colnames(back_inloop2)[1]),"background",sep = "_"),
              paste(gsub('.{2}$','',colnames(back_inloop2)[1]),"background_others",sep = "_"),
              paste(gsub('.{2}$','',colnames(bait_inloop2)[1]),"bait",sep = "_"),
              paste(gsub('.{2}$','',colnames(bait_inloop2)[1]),"bait_others",sep = "_")      
            )
            
            
            
            #===============================================#
            #---------Start Fisher Loop   num1 = 0 ---------#
            #===============================================#
            
            
            data3 <- in_df[,6:9]
            
            
            k_row <- nrow(data3)
            
            
            Fisher_P <-data.frame()
            
            
            for (k in 1:k_row) {
              
              #~~~#
              tryCatch({   
                #~~~#
                
                
                datainput <- matrix(as.numeric(c(data3[k,1],data3[k,2],data3[k,3],data3[k,4])), 2,
                                    dimnames = list(c("protein","others"),
                                                    c("background","bait")
                                    )
                )
                
                dataoutput <- fisher.test(datainput, alternative = "less")
                
                
                Fisher_P[k,1] <- dataoutput$p.value
                
                #~~~#  
              }, error=function(e){print(paste("=======","Error occurred in",result_out[i],data2[k,1],"=======",sep = " "))})  
              #~~~# 
            }
            
            colnames(Fisher_P) <- "Fisher_P"
            
            result <- cbind(in_df,Fisher_P)
            
            result_cutoff <- result[which(result$Fisher_P<input$pvalue),]
            
            
            
            #collect protein ID after cutoff
            
            
            medi1 <- data.frame(result_cutoff[,2])
            
            colnames(medi1) <- "V1"
            
            medi1 <- medi1 %>% filter(!grepl('\\*', V1))
            
            result_cutoff_protein_list[i] <- list(medi1[,1])
            
            names(result_cutoff_protein_list)[i] <- paste(gsub('.{2}$','',colnames(bait_inloop2)[1]), "_to_", gsub('.{2}$','',colnames(back_inloop2)[1]) ,sep = "")
            
          }
          
          
          
          result_cutoff_protein_list
          
          
        }
      
      
      
    }
    
    
  })
  
  #generate bait target name, pass to Upset plot output function
  getBait_target_name <- reactive({
    
    if(input$goButton!= 0){
      
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "csv", "Please upload a csv file"))
      
      bait_target_name1 <- gsub(".csv","",file)
      
      bait_target_name1
    }
  })
  
  #list of prey ---> chart
  ListofPrey.table  <- reactive({
    
    if(input$goButton!= 0 && input$num1 == 0){
      
      list1 <- getData1()
      
      dat <- data.table::rbindlist(
        lapply(list1, function(x) data.table::data.table(t(x))),fill = TRUE) %>% t() %>% 
        data.frame(row.names = seq(1:max(lengths(list1)) )) %>% 
        purrr::set_names(names(list1))
      
      dat
    }
  })
  
  
  #================================== Result show in interface =================================#
  
  #-----------------------------------List of Prey / fisher table ------------------------------------#
  which.table  <- reactive({
    
    if(input$goButton!= 0 && input$num1 == 0){
      
      show.table <- ListofPrey.table()
      
      show.table
      
    } else if (input$goButton!= 0 && input$num1 != 0) {
      
      show.table <- getData1()
      
      #digit number
      show.table <- data.frame(format(show.table, digit = 9))
      
      show.table
      
    }
    
    
  })
  
  output$contents  <-  renderTable({
    
    data <- which.table()
    
    data
    
  })
  
  
  
  #----------------------------------- upset plot reactive ------------------------------------#
  Get_upsetplot <-  reactive({
    
    
      
      list1 = getData1()
      
      bait_target_name2 <- getBait_target_name()
     
      p1 <- upset(fromList(list1),  
                  
                  nsets = 30,     
                  
                  nintersects = 40, 
                  
                  order.by = "freq", 
                  
                  keep.order = TRUE,
                  
                  mb.ratio = c(0.45,0.55),   
                  
                  text.scale = 1.9,  
                  
                  #mainbar.y.label = paste("intersections ", "(", bait_target_name2, ")", sep = "")   #error/warning occur after R 4.1
      )
      
      p1
    
     
  })
  
  #----------------------------------- show upset plot ------------------------------------#
  
  
  output$contents2 <- renderPlot({
    
    if(input$goButton!= 0 && input$num1 == 0){
    
    pic <- Get_upsetplot()
    
    pic 
    
    } else {print("Please set '0' for parameter #2") 
      
      x<-matrix(1,1)
      p2 <- plot(x,ann = F, bty = "n", xaxt = "n", yaxt ="n",col ="white")
      title("Please set '0' for parameter #2 for Upset plotting after clicking ‘Run’ ")
      
    }  
      
      
    
  })
  
  
  
  
  
  #-------------------------- interaction chart -------------------------#
  
  
  inter_chart  <-  reactive({
    
    if(input$goButton!= 0 && input$num1 == 0){
      
      list1 = getData1()
      
      df_inter <- get.venn.partitions(list1)
      
      for (i in 1:nrow(df_inter)) df_inter[i,'elements'] <- paste(df_inter[[i,'..values..']], collapse = ', ')
      
      nn <- ncol(df_inter)
      
      df__inter_fin <- df_inter[-c(nn-3,nn-2)]
      
      
      
      df__inter_fin
      
    } 
    
    
    
  })
  
  output$contents3  <-  renderTable({
    
    data_int = inter_chart()
    
    data_int
    
  }, width = "400px")
  
  
  #------------------------- download function ------------------------#
  
  output$tableDown <- downloadHandler(
    
    filename = paste("result_Fisher",Sys.Date(),".txt",sep = ""),
    
    content = function(file){
      
      vroom_write(which.table(), file, delim = "\t")
      
    }
  )
  
  
  output$upsetplotDown <- downloadHandler( 
    
    filename = paste("upset_plot",Sys.Date(),".pdf",sep = ""),
    
    content = function(file) {
      pdf(file ,height = 16, width = 22)
      print(Get_upsetplot())
      dev.off()                                      
      
      
    }
  )
  
  
  output$tableDown_intersect <- downloadHandler(
    
    filename = paste("intersection_chart",Sys.Date(),".txt",sep = ""),
    
    content = function(file){
      
      vroom_write(inter_chart(), file, delim = "\t")
      
    }
  )
  
  
}


####################
shinyApp(ui, server)