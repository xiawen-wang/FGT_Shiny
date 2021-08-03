#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/

#######################################################################
# This is a Shiny script for an interactive volcano plot based on the #
# statistic analysis results of Limma. The differentially expressed   #
# genes can be easily identify using a volcano plot. Here, genes      #
# upregulated or downregulated are highlighted by different colours   #
# interactive options are provided for significance cutoffs, so the   #
# user can find the genes of interested.  Besides, Font sizes and     #
# point sizes can also be changed for better visual experiences.      #
# Moreover, user can click on the points of interest or select        #
# interested data area to display more detailed point information.    #
# Finally, a table of all data information is provided in a seperate  #
# panel for easy data access. Genes can be easily sorted or searched. #
#                                                                     #
# Note: this warning message might be seen when the app is run:       #
# "Warning: Removed 45094 rows containing missing values              #
# (geom_text_repel)". This is because only upregulated and            #
# downregulated gene names are displayed for the plot. This does not  #  
# affect the functionality or accuracy of the plot.                   #
#######################################################################

#### Reference for the volcano plot:
#    https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

#### References for the added functions:
#    Shiny skeleton:http://shiny.rstudio.com/
#    Selection/Clicking points: https://shiny.rstudio.com/gallery/plot-interaction-selecting-points.html
#    Display table for clicks or selections: https://shiny.rstudio.com/articles/render-table.html
#    Navigation bar: https://shiny.rstudio.com/gallery/navbar-example.html
#    Data table: https://rstudio.github.io/DT/shiny.html

#######################################################################


# The following packages are needed for the script to run
library(shiny)
library(ggplot2)
library(ggrepel)
library(plyr)
library(markdown)
library(DT)

#load the example data
load("volcano.Rdata")


## Volcano plot
# make a data frame with all the data needed for ploting
vol.data=data.frame(Gene=myresults$Symbol,logFC=myresults$logFC,P.Value=myresults$P.Value,adj.P.Val=myresults$adj.P.Val,B=myresults$B)

# make a new column for -log10(adj.p.Val), which will be ploted on the y axis
vol.data$log10fdr <- NA
vol.data$log10fdr <- -log10(vol.data$adj.P.Val)


# Define UI for application that draws a Volcano plot
ui <- fluidPage(

    # CSS for text style
    tags$head(
        tags$style("
                                     
             .navbar{
             font-size: 16px;
             font-weight: bold;
             }
             
             #upreg{
             color: #d42020;
             font-size: 15px;
             font-weight: bold;
             }
             
             #downreg{
             color: #2080d4;
             font-size: 15px;
             font-weight: bold;
             }"
        )
    ),
    
    # Navigation bar
    
    navbarPage("Differential Expression Analysis",
                tabPanel("Volcano Plot",
                    
                    # Control elements for point selection     
                    fluidRow(
                        column(3,
                        sliderInput("log10fdr",
                                   "-log10(FDR adj. P-value):",
                                   min = 0,
                                   max = 3,
                                   value = 1.3,
                                   step = 0.05)
                        ),
                        column(3,
                        sliderInput("log2fc",
                                   "Log2FC:",
                                   min = 0,
                                   max = 4,
                                   value = 1,
                                   step = 0.1)
                        ),
                        column(3,
                        sliderInput("fontsize",
                                   "Font size:",
                                   min = 8,
                                   max = 24,
                                   value = 14,
                                   step = 1)
                        ),
                        column(3,
                        sliderInput("pointsize",
                                   "Point size:",
                                   min = 1,
                                   max = 10,
                                   value = 4,
                                   step = 1)
                        )
                    ),
                    hr(),
                    
                    fluidRow(
                    
                    # The plot
                        column(7,   
                            plotOutput("plot1",width="100%", 
                                       height="600", 
                                       click = "plot1_click",
                                       brush = brushOpts(id = "plot1_brush")
                            ),
                        ),
                    
                    # column for text display and selected data table display 
                        column(3,
                            textOutput("upreg"),
                            textOutput("downreg"),
                            
                            hr(),
                            
                            h4("Points near click"),
                            tableOutput("click_info"),
                        
                        
                            h4("Points within selection"),
                            tableOutput("brush_info")
                        )
                    )
                ),
               
               # a separate panel to display all data
                tabPanel("Complete Data Table",
                    DT::dataTableOutput("table")
                )
            )
)



server <- function(input, output) {

    
    output$plot1 <- renderPlot({

        # make a column to differentiate upregulated, downregulated and no change genes 
        vol.data$diffexpressed <- "NO"
        vol.data$diffexpressed[vol.data$logFC > input$log2fc & (vol.data$log10fdr) > input$log10fdr]<-"UP"
        vol.data$diffexpressed[vol.data$logFC < (-input$log2fc) & (vol.data$log10fdr) > input$log10fdr]<-"DOWN"
        
        # make a colum to have gene names only for upregulated or downregulated genes
        vol.data$glabel <- NA
        vol.data$glabel[vol.data$diffexpressed!="NO"]<-vol.data$Gene[vol.data$diffexpressed!="NO"]
        
        # some information about the cutoffs
        plottitle <- paste("P-Value cutoff: ", (10^(-(input$log10fdr))), "\nLog2 fold change cutoff: +", input$log2fc,"& -", input$log2fc)
        
        # plot the volcano plot
        ggplot(vol.data,aes(x=logFC, y=log10fdr,col=diffexpressed, label=glabel)) +
        theme_minimal() + 
        ggtitle(plottitle) +
        geom_point(size=input$pointsize, alpha=0.25, stroke = 1) + 
        
        # gene name display do not overlap
        geom_text_repel(size=(input$fontsize)*0.352777778) +
        
        # axis names and breaks
        scale_x_continuous(name="Fold change (log2)", breaks = seq(-4, 4, by = 1)) +
        scale_y_continuous(name="-log10(FDR-adj. P-value") +
        
        # horizontal line for the p value cutoff, vertical line for the fold change cutoff
        geom_vline(xintercept=c((input$log2fc)*-1, input$log2fc), col="black", linetype=2) +
        geom_hline(yintercept=(input$log10fdr), col="red", linetype=5)+
        geom_vline(xintercept = c(-4,4), linetype = 2, alpha=0)+
        
        # color scheme for the points
        scale_color_manual(name="Differential expression", values=c("NO" = "#999999", "UP" = "#d42020", "DOWN" = "#2080d4")) +
        
        # font sizes
        theme(plot.title = element_text(size=input$fontsize), 
              legend.position="bottom", 
              axis.text.x = element_text(size=input$fontsize), 
              axis.title.x = element_text(size=input$fontsize), 
              axis.text.y = element_text(size=input$fontsize), 
              axis.title.y = element_text(size=input$fontsize),
              legend.text=element_text(size=input$fontsize)) 
       
        
    })
    
    # display table of approximate point information on click
    output$click_info <- renderTable({
        # Because it's a ggplot2, xvar or yvar are not needed; if this
        # were a base graphics plot, we'd need those.
        nearPoints(vol.data, input$plot1_click, addDist = FALSE)
    })
    
    # display point information in a selected area
    output$brush_info <- renderTable({
        brushedPoints(vol.data, input$plot1_brush)
    })
    
    # all data information for the second penal
    output$table <- DT::renderDataTable({
        DT::datatable(myresults)
    })
    
    # display text information for upregulated gene number
    output$upreg <- renderText({
        vol.data$diffexpressed <- "NO"
        vol.data$diffexpressed[vol.data$logFC > input$log2fc & (vol.data$log10fdr) > input$log10fdr]<-"UP"
        vol.data$diffexpressed[vol.data$logFC < (-input$log2fc) & (vol.data$log10fdr) > input$log10fdr]<-"DOWN"
        upreg <- paste("Upregulated genes: ", count(vol.data$diffexpressed=="UP")[2,2])
        upreg
        
    })
    
    # display text information for downregulated gene number
    output$downreg <- renderText({
        vol.data$diffexpressed <- "NO"
        vol.data$diffexpressed[vol.data$logFC > input$log2fc & (vol.data$log10fdr) > input$log10fdr]<-"UP"
        vol.data$diffexpressed[vol.data$logFC < (-input$log2fc) & (vol.data$log10fdr) > input$log10fdr]<-"DOWN"
        downreg <- paste("Downregulated genes: ", count(vol.data$diffexpressed=="DOWN")[2,2])
        downreg
        
    })
    

}

# Run the application 
shinyApp(ui = ui, server = server)
