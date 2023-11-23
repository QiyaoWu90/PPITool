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

# library(dplyr)
# library(UpSetR)
# library(VennDiagram)



###################################################################################
#================================     UI      ====================================#
###################################################################################

ui <- fluidPage(
  
 
  #navbar icon pic
  navbarPage(title=div(img(src="MyIcon.png",height = 25, width = 30), "PPI Tool"),
          
            
             
             #background pic
             tags$style(
               'div[data-value="Home"]{ 
                 width: 100%; height: 950px; 
                 background-image: url("background_img_mainpage.jpg");
               }'
             ),
             
             #============================== Main page ================================#  
             tabPanel("Home",
                      
                      
                      titlePanel(h1("Welcome to PPI Tool",align="center",style="color: skyblue")), 
                      
                      titlePanel(h4("PPI Tool uses a probabilistic method based on fisher' exact test ",align="center")),   
                      titlePanel(h4("screens out potential significants protein-protein interactors to your bait of interests",align="center")),    
                      titlePanel(h4("It demonstrates better performance and higher true positive",align="center")), 
                      titlePanel(h4("gives you more potential choices for proteomic interaction studies",align="center")), 
                      titlePanel(h3("For the usage of PPI Tool, check 'Documentation' tab",align="center",style="color: blue")),
                      titlePanel(h3("To start the prediction, click 'analysis' tab",align="center",style="color: blue")),
                   
                       
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
                        img(src="FAFUicon.jpg", height = 100, width = 100),
                        br(),
                        br(),
                        
                        p(""),
                        br(),
                        p("The PPI Tool is created by Qiyao Wu and Xueyang Zhang Ph.D"),
                        br(),
                        
                        
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
      
      #提取数据列，前五列为char信息，从6开始为numeric data
      data_num <- data.frame(data1[,c(6:col_number)])
      
      #提取bait，数据前c(1:input$repeatsn)为bait  (input$repeats = 重复数) 
      bait <- data_num[,c(1:input$repeats)]
      
      #提取背景  (5是数据前char 信息列的列数)
      back <- data_num[,c((input$repeats +1 ):(col_number-5))]
      
      #多少组背景（每个背景有input$repeats个重复）
      n <- ncol(back)/input$repeats
      
      #空list用于收集不同背景的比对结果
      result_cutoff_protein_list <- list()
      
      
      #===================================================================================================#
      #num1 为非0时，计算单个bait vs background并显示结果
      
      if(input$num1 != 0) {
        
        i <- input$num1  
        
        
        #根据repeats 读取每组对比中的background数据集
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
        
        #提取char信息
        in_df <- df[,c(1:5)]
        
        #背景总和(所有重复，暂定)
        background_total <- sum(back_inloop)
        
        #bait总和(所有重复，暂定)
        bait_total <- sum(bait)
        
        
        #backg. 和bait总数比大小，决定情况
        if (background_total < bait_total)
          
        {
          back_inloop2 <- back_inloop
          bait_inloop2 <- bait
          
          
          
          background_total_loop2 <- background_total
          bait_total_loop2       <- bait_total
          
        } 
        
        #如果bait比background小，要impute所有0为1
        if (background_total > bait_total)
          
        {
          
          df2 <- cbind(in_df,bait)
          df2_original <- df2      #用来比对impute前后
          
          row_num <- nrow(df2)
          
          rows_of_zeros <- 0
          
          #计数全是0的row
          for (x in 1:row_num) 
            
          {
            matrix1 <- as.matrix(df2[x,c(6:8)])
            
            
            
            if (matrix1[1,1] ==0 &  matrix1[1,2] ==0 & matrix1[1,3] ==0 ) {
              rows_of_zeros = rows_of_zeros +1
            }
            
          }
          
          #计算bait 和background 的差,与全是0的row相差几倍
          fold1 <- (background_total - bait_total)/rows_of_zeros
          
          #定义impute数,给全是0的row impute一个数,来使bait总和超过background （最多不超过4）
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
          
          #根据上述结果对df2进行impute
          for (x in 1:row_num) 
          {
            matrix1 <- as.matrix(df2[x,c(6:8)])
            
            if (matrix1[1,1] ==0 &  matrix1[1,2] ==0 & matrix1[1,3] ==0 )
              
            { 
              matrix1[1,1] <- impute_num      #impute number
              df2[x,c(6:8)]    <-    matrix1  
              df2[x,2] <- paste(df2[x,2],sep = "")   #在majority protein上加星号，值是impute_star,(先不加)，注意majority的位置是[x,2]
            }
          }
          
          #重新提取char信息（更新带有impute星号的majority protein name）  
          in_df <- df2[,c(1:5)] 
          
          back_inloop2 <- back_inloop
          bait_inloop2 <- df2[,c(6:8)]
          
          
          
          background_total_loop2 <- background_total
          bait_total_loop2       <- sum(df2[,c(6:8)]) 
          
        }
        
        
        
        
        #生成新col信息
        
        #1 添加各个row background 数据（暂定重复总和）
        in_df$background <- apply(back_inloop2,1,sum)
        
        
        #2添加background othesr列
        in_df$background_others <- background_total_loop2 - in_df$background
        
        #bait row总和  
        in_df$baited <-  apply(bait_inloop2,1,sum)
        
        #3计算bait 
        in_df$bait <- in_df$baited -  in_df$background
        
        #去掉baited列
        in_df <- subset(in_df,select = -c(baited))
        
        #负值变为0
        in_df$bait[in_df$bait < 0 ] = 0
        
        
        #4计算bait_others 
        in_df$bait_others <- bait_total_loop2 - background_total_loop2 - in_df$bait
        
        #负值变为0
        in_df$bait_others[in_df$bait_others < 0 ] = 0
        
        
        
        #重命名col
        colnames(in_df)[6:9] <- c(
          paste(gsub('.{2}$','',colnames(back_inloop2)[1]),"background",sep = "_"),
          paste(gsub('.{2}$','',colnames(back_inloop2)[1]),"background_others",sep = "_"),
          paste(gsub('.{2}$','',colnames(bait_inloop2)[1]),"bait",sep = "_"),
          paste(gsub('.{2}$','',colnames(bait_inloop2)[1]),"bait_others",sep = "_")      
        )
        
        
        #================================#
        #---------开始fisher循环---------#
        #================================#
        
        
        data3 <- in_df[,6:9]
        
        #读取行数
        k_row <- nrow(data3)
        
        #空数据集收集P值
        Fisher_P <-data.frame()
        
        
        for (k in 1:k_row) {
          
          #~~~#
          tryCatch({   #用tryCatch功能夹住方程，如果因无法算P产生错误则跳过，生成NA
            #~~~#
            
            #构建fisher矩阵
            datainput <- matrix(as.numeric(c(data3[k,1],data3[k,2],data3[k,3],data3[k,4])), 2,
                                dimnames = list(c("protein","others"),
                                                c("background","bait")
                                )
            )
            #Fisher Test
            dataoutput <- fisher.test(datainput, alternative = "less")
            
            #收集P值
            Fisher_P[k,1] <- dataoutput$p.value
            
            #~~~#  #console内提示错误
          }, error=function(e){print(paste("=======","Error occurred in",result_out[i],data2[k,1],"=======",sep = " "))})  #用tryCatch功能夹住方程，如果因无法算P产生错误则跳过，生成NA-尾部
          #~~~# 
        }
        
        colnames(Fisher_P) <- "Fisher_P"
        
        result_cutoff_protein_list <- data.frame(cbind(in_df,Fisher_P))
        
        
        
        result_cutoff_protein_list
        
        
      } else
        if (input$num1 == 0) {
          
          #===================================================================================================================#
          #当num1 = 0 则开始构建不同bait-background的df
          for (i in 1:n) {
            
            
            #根据repeats 读取每组对比中的background数据集
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
            
            #提取char信息
            in_df <- df[,c(1:5)]
            
            #背景总和(所有重复，暂定)
            background_total <- sum(back_inloop)
            
            #bait总和(所有重复，暂定)
            bait_total <- sum(bait)
            
            
            #backg. 和bait总数比大小，决定情况
            if (background_total < bait_total)
              
            {
              back_inloop2 <- back_inloop
              bait_inloop2 <- bait
              
              
              
              background_total_loop2 <- background_total
              bait_total_loop2       <- bait_total
              
            } 
            
            #如果bait比background小，要impute所有0为1
            if (background_total > bait_total)
              
            {
              
              df2 <- cbind(in_df,bait)
              df2_original <- df2      #用来比对impute前后
              
              row_num <- nrow(df2)
              
              rows_of_zeros <- 0
              
              #计数全是0的row
              for (x in 1:row_num) 
                
              {
                matrix1 <- as.matrix(df2[x,c(6:8)])
                
                
                
                if (matrix1[1,1] ==0 &  matrix1[1,2] ==0 & matrix1[1,3] ==0 ) {
                  rows_of_zeros = rows_of_zeros +1
                }
                
              }
              
              #计算bait 和background 的差,与全是0的row相差几倍
              fold1 <- (background_total - bait_total)/rows_of_zeros
              
              #定义impute数,给全是0的row impute一个数,来使bait总和超过background （最多不超过4）
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
              
              #根据上述结果对df2进行impute
              for (x in 1:row_num) 
              {
                matrix1 <- as.matrix(df2[x,c(6:8)])
                
                if (matrix1[1,1] ==0 &  matrix1[1,2] ==0 & matrix1[1,3] ==0 )
                  
                { 
                  matrix1[1,1] <- impute_num      #impute number
                  df2[x,c(6:8)]    <-    matrix1  
                  df2[x,2] <- paste(df2[x,2],sep = "")   #在majority protein上加星号，值是impute_star,(先不加)，注意majority的位置是[x,2]
                }
              }
              
              #重新提取char信息（更新带有impute星号的majority protein name）  
              in_df <- df2[,c(1:5)] 
              
              back_inloop2 <- back_inloop
              bait_inloop2 <- df2[,c(6:8)]
              
              
              
              background_total_loop2 <- background_total
              bait_total_loop2       <- sum(df2[,c(6:8)]) 
              
            }
            
            
            
            
            #生成新col信息
            
            #1 添加各个row background 数据（暂定重复总和）
            in_df$background <- apply(back_inloop2,1,sum)
            
            
            #2添加background othesr列
            in_df$background_others <- background_total_loop2 - in_df$background
            
            #bait row总和  
            in_df$baited <-  apply(bait_inloop2,1,sum)
            
            #3计算bait 
            in_df$bait <- in_df$baited -  in_df$background
            
            #去掉baited列
            in_df <- subset(in_df,select = -c(baited))
            
            #负值变为0
            in_df$bait[in_df$bait < 0 ] = 0
            
            
            #4计算bait_others 
            in_df$bait_others <- bait_total_loop2 - background_total_loop2 - in_df$bait
            
            #负值变为0
            in_df$bait_others[in_df$bait_others < 0 ] = 0
            
            
            
            #重命名col
            colnames(in_df)[6:9] <- c(
              paste(gsub('.{2}$','',colnames(back_inloop2)[1]),"background",sep = "_"),
              paste(gsub('.{2}$','',colnames(back_inloop2)[1]),"background_others",sep = "_"),
              paste(gsub('.{2}$','',colnames(bait_inloop2)[1]),"bait",sep = "_"),
              paste(gsub('.{2}$','',colnames(bait_inloop2)[1]),"bait_others",sep = "_")      
            )
            
            
            
            #================================#
            #---------开始fisher循环---------#
            #================================#
            
            
            data3 <- in_df[,6:9]
            
            #读取行数
            k_row <- nrow(data3)
            
            #空数据集收集P值
            Fisher_P <-data.frame()
            
            
            for (k in 1:k_row) {
              
              #~~~#
              tryCatch({   #用tryCatch功能夹住方程，如果因无法算P产生错误则跳过，生成NA
                #~~~#
                
                #构建fisher矩阵
                datainput <- matrix(as.numeric(c(data3[k,1],data3[k,2],data3[k,3],data3[k,4])), 2,
                                    dimnames = list(c("protein","others"),
                                                    c("background","bait")
                                    )
                )
                #Fisher Test
                dataoutput <- fisher.test(datainput, alternative = "less")
                
                #收集P值
                Fisher_P[k,1] <- dataoutput$p.value
                
                #~~~#  #console内提示错误
              }, error=function(e){print(paste("=======","Error occurred in",result_out[i],data2[k,1],"=======",sep = " "))})  #用tryCatch功能夹住方程，如果因无法算P产生错误则跳过，生成NA-尾部
              #~~~# 
            }
            
            colnames(Fisher_P) <- "Fisher_P"
            
            result <- cbind(in_df,Fisher_P)
            
            result_cutoff <- result[which(result$Fisher_P<input$pvalue),]
            
            #收集cutoff后的protein ID 到总表内
            result_cutoff_protein_list[i] <- list(result_cutoff[,2])
            
            names(result_cutoff_protein_list)[i] <- paste(gsub('.{2}$','',colnames(bait_inloop2)[1]), "_to_", gsub('.{2}$','',colnames(back_inloop2)[1]) ,sep = "")
            
          }
          
          
          
          result_cutoff_protein_list
          
          
        }
      
      
      
    }
    
    
  })
  
  #生成bait target name 传导给Upset plot output功能
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
  
  #把list of prey变成表格
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
     
      p1 <- upset(fromList(list1),  #转换兼容数据形式
                  
                  nsets = 30,     # 绘制的最大集合个数 - 对应表格col，图下-柱形图+点图
                  
                  nintersects = 40, #绘制的最大交集个数，NA则全部绘制 - 对应表格rol，图上-柱形图
                  
                  order.by = "degree", # 矩阵中的交点是如何排列的。 "freq"根据交集个数排序，"degree" 根据基因交集密集度
                  
                  keep.order = TRUE,
                  
                  mb.ratio = c(0.45,0.55),   # 左侧和上方条形图的比例关系
                  
                  text.scale = 1.9, # 文字标签的大小
                  
                  mainbar.y.label = paste("intersections ", "(", bait_target_name2, ")", sep = "")
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
  
  
  
  
  
  #------------------------- interaction chart ------------------------#
  
  
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



shinyApp(ui, server)