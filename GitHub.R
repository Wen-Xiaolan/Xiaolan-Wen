
#=========================================cancer number==============================================
library(dplyr)
PAF<- read.xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_age.xlsx")
PAF<-PAF[,c(5,6)]
PAF <- PAF %>%
  distinct()

library(readr)
case <- read_csv("~/Desktop/广东省肿瘤数据/cancer_all.csv")
case<-case[,-c(1,2,8)]


case <- case %>%
  filter(age != "f1" & age != "f2" & age != "f3" & age != "f4") 

case$ICD10[case$ICD10 %in% c("C18","C19-C20")]<-"C19-C20"

case<-case%>%
  left_join(PAF,by="ICD10")

#case <- case %>%
#filter(ICD10 == "ALL")

case<-aggregate(cbind(inc_all,mort_all)~sex+kind+ICD10,case,sum)

case$inc_all <- as.integer(case$inc_all)
case$mort_all <- as.integer(case$mort_all)


case <- case %>%
  filter(kind == 0)
case <- case %>%
  filter(sex == 0)

df_sorted <- case %>%
  arrange(desc(mort_all))

data<- read.xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_case_all_exposure_ICD10(all_age).xlsx")
data <- data %>%
  filter(kind == 2)
data <- data %>%
  filter(sex == 0)

df_sorted <- data %>%
  arrange(desc(mort_attribute))

#====================================================meta分析==========================================
install.packages("readxl")
install.packages("meta")
library(readxl)
library(meta)

data <- read_excel("~/Desktop/广东省肿瘤数据/暴露/感染因素-暴露.xlsx", sheet = "EB-NPC")

data$exposure_rate <- data$case / data$total

data$SE <- sqrt((data$exposure_rate * (1 - data$exposure_rate)) / data$total)

#data <- data[order(data$Study), ]

meta_analysis <- metaprop(event = data$case, n = data$total, data = data, sm = "PLO", method = "Inverse")

summary(meta_analysis)

# 确保将 Study 列赋值给 meta_analysis 对象的 studlab 属性
meta_analysis$studlab <- data$Study  # 将 Study 列作为研究标签

tiff("~/Desktop/广东省肿瘤数据/暴露/EBV_nanolaryngeal-forest_plot.tiff", width = 12, height = 18, units = "in", res = 1000)

# 绘制森林图
par(family = "Times New Roman")
forest(meta_analysis, 
       xlab = "Prevalence",                  # x轴标签
       slab = meta_analysis$studlab,       # 使用 Study 列作为研究名称
       cex = 0.7,                          # 字体大小
       alim = c(0, 1),                     # 设置x轴范围
       refline = 0.5,                      # 添加参考线
       showweights = FALSE,                # 不显示加权数据
       showaxis = TRUE,                    # 显示轴
       xlab.cex = 0.3                         # 调整x轴标签字体大小
)


dev.off()

#==========================平均风险方法计算饮食因素的PAF==========================================

install.packages("readxl")
install.packages("boot")

library(readxl)
library(dplyr)
library(writexl)
library(boot)
library(openxlsx)

RR <- read_excel("~/Desktop/广东省肿瘤数据/暴露/广东暴露率及RR（sex）.xlsx", sheet = "RR-饮食")                                                                        
exposure <- read_excel("~/Desktop/广东省肿瘤数据/暴露/广东暴露率及RR（sex）.xlsx", sheet = "暴露率-饮食")

RR$se<-paste(RR$sex,RR$exposure,RR$age,RR$kind)
exposure$se<-paste(exposure$sex,exposure$exposure,exposure$age,exposure$kind)

exposure<-exposure[,-c(1:4)]

merged_data <- full_join(RR, exposure, by = "se")

data<-merged_data

data<-data[,-c(10,11)]

# 找出包含NA的行
rows_with_na <- data[!complete.cases(data), ]
#print(rows_with_na)

# 1. 计算风险 (risk)
data$risk <- exp(log(data$RR) * data$rate)

# 2. 计算PAF
data$PAF <- (data$risk - 1) / data$risk

# 3. 计算Delta法方差 (var_PAF)
# 假设 var_P（暴露率的方差）为零，暴露率P在此为rate列
# var_beta = (RR_upper - RR_lower) / (2 * RR) 这是一个近似值，用于估算RR的标准误差。

# 计算RR的标准误差 (这里使用对数变换的标准误差)
data$var_beta <- ((data$RR_upper - data$RR_lower) / (2 * data$RR))^2

# 计算PAF的方差 (var_PAF)
data$var_PAF <- ((data$risk - 1)^2 * 0 + (data$rate * data$risk)^2 * data$var_beta) / ((data$rate * (data$risk - 1) + 1)^4)

# 4. 计算PAF的标准误差 (SE_PAF)
data$se_PAF <- sqrt(data$var_PAF)

# 5. 计算95%置信区间
data$PAF_lower <- data$PAF - 1.96 * data$se_PAF
data$PAF_upper <- data$PAF + 1.96 * data$se_PAF


write.xlsx(data, "~/Desktop/广东省肿瘤数据/PAF2019/attribute_diatery.xlsx")

#======================================Levin法的PAF==========================================
install.packages("readxl") 
install.packages("dplyr")  
install.packages("writexl") 

library(readxl)
library(dplyr)
library(writexl)
library(openxlsx)

RR <- read_excel("~/Desktop/数据/广东省肿瘤归因/暴露/广东暴露率及RR（sex）.xlsx", sheet = "RR")                                                                        
exposure <- read_excel("~/Desktop/数据/广东省肿瘤归因/暴露/广东暴露率及RR（sex）.xlsx", sheet = "暴露率-sex+kind")

RR$se<-paste(RR$sex,RR$exposure,RR$age,RR$kind)
exposure$se<-paste(exposure$sex,exposure$exposure,exposure$age,exposure$kind)

exposure<-exposure[,-c(1:4)]

merged_data <- full_join(RR, exposure, by = "se")

data<-merged_data

# 找出包含NA的行
rows_with_na <- data[!complete.cases(data), ]
# print(rows_with_na)

# 删除包含NA的行
# data_clean <- data[complete.cases(data), ]
# data<-data_clean


# 计算PAF及其95%置信区间的函数
calculate_PAF_and_CI <- function(RR, RR_lower, RR_upper, P) {
  # 计算beta (ln(RR))
  beta <- log(RR)
  beta_lower <- log(RR_lower)
  beta_upper <- log(RR_upper)
  
  # 计算beta的方差（Var(beta))
  beta_se <- (beta_upper - beta_lower) / (2 * 1.96) 
  var_beta <- beta_se^2
  exp_beta <- exp(beta)
  PAF <- ((exp_beta - 1) * P) / (P * (exp_beta - 1) + 1)
  var_P <- 0  # 假设暴露率P的方差为零（如果P是整个群体的暴露率）
  var_PAF <- ((exp_beta - 1)^2 * var_P + (P * exp_beta)^2 * var_beta) / ((P * (exp_beta - 1) + 1)^4)
  se_PAF <- sqrt(var_PAF)
  PAF_lower <- PAF - 1.96 * se_PAF
  PAF_upper <- PAF + 1.96 * se_PAF
  
  return(list(PAF = PAF, PAF_lower = PAF_lower, PAF_upper = PAF_upper, SE_PAF = se_PAF))
}


# 使用apply函数来批量计算每行的PAF和95%置信区间
results <- apply(data, 1, function(row) {
  # 确保数据是数值型
  RR <- as.numeric(row["RR"])
  RR_lower <- as.numeric(row["RR_lower"])
  RR_upper <- as.numeric(row["RR_upper"])
  rate <- as.numeric(row["rate"])
  
  # 检查是否有任何值是NA或无效的
  if (is.na(RR) || is.na(RR_lower) || is.na(RR_upper) || is.na(rate)) {
    return(c(NA, NA, NA, NA))  # 如果有NA，返回NA
  }
  
  # 调用PAF计算函数
  result <- calculate_PAF_and_CI(RR, RR_lower, RR_upper, rate)
  
  # 返回计算结果
  return(c(result$PAF, result$PAF_lower, result$PAF_upper, result$SE_PAF))
})

# 转换计算结果为数据框
results_df <- as.data.frame(t(results))
colnames(results_df) <- c("PAF", "PAF_lower", "PAF_upper", "SE_PAF")

# 将计算结果合并到原始数据框
data <- cbind(data, results_df)

# 要保留4位小数的列名
# cols_to_round <- c("PAF", "PAF_lower", "PAF_upper", "SE_PAF")  
# 使用round函数保留4位小数
# data[cols_to_round] <- round(data[cols_to_round], 4)

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/PAF2019/attribute_behaviour+metabolic.xlsx")

#========================糖尿病和BMI的PAF==========================================

data$kasi <- paste(data$kind,data$age,data$sex,data$ICD10)

# 找出每个 age 中 BMI 和 diabetes 的 PAF 较小的那一行
to_remove <- data %>%
  filter(exposure %in% c("Excess bodyweight", "Diabetes")) %>%
  group_by(kasi) %>%
  slice_min(order_by = PAF, n = 1, with_ties = FALSE) %>%
  ungroup()

# 从原数据中移除这些行
result <- anti_join(data, to_remove, by = c("kasi", "exposure", "PAF"))

data <- result

data <- data %>%
  mutate(exposure = if_else(exposure == "Excess bodyweight", "Diabetes", exposure))

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_behaviour+metabolic1.xlsx")

#========================合并膳食、感染和其他因素的PAF==========================================

data1 <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_behaviour+metabolic1.xlsx")
data2 <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_diatery.xlsx")
data3 <- read_excel("~/Desktop/数据/广东省肿瘤归因/暴露/广东暴露率及RR（sex）.xlsx", sheet = "暴露因素-患者中的率") 


data1<-data1[,c(1:6,12:14)]
data2<-data2[,c(1:6,15,19:20)]


#colnames(data2)[colnames(data2) == "se_PAF"] <- "SE_PAF"

data<-rbind(data1,data2,data3)

# 将 PAF 列中负数的值改为 0
data$PAF[data$PAF < 0] <- 0
data$PAF_lower[data$PAF_lower < 0] <- 0


write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_age.xlsx")

#========================风险对的PAF及归因病例（分年龄、分城乡、分性别）==========================================

library(readr)

PAF<- read.xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_age.xlsx")


case <- read_csv("~/Desktop/数据/广东省肿瘤归因/分析结果/cancer_all.csv")
case<-case[,-c(1,2,8)]

case$age[case$age %in% c("f5","f6")]<-"20-29"
case$age[case$age %in% c("f7","f8")]<-"30-39"
case$age[case$age %in% c("f9","f10")]<-"40-49"
case$age[case$age %in% c("f11","f12")]<-"50-59"
case$age[case$age %in% c("f13","f14")]<-"60-69"
case$age[case$age %in% c("f15","f16","f17","f18")]<-"70-"

case <- case %>%
  filter(age != "f1" & age != "f2" & age != "f3" & age != "f4") 

case1<-aggregate(cbind(inc_all,mort_all)~ICD10+sex+kind+age,case,sum)

case<-case1

case$ICD10[case$ICD10 %in% c("C00","C01-C02","C03-C06","C07-C08","C09","C10","C12-C13","C14")]<-"C00-C14,except C11"
case$ICD10[case$ICD10 %in% c("C51","C52")]<-"C51-C52"
case1<-aggregate(cbind(inc_all,mort_all)~ICD10+sex+kind+age,case,sum)
case<-case1

case$iska<-paste(case$ICD10,case$sex,case$kind,case$age)
case<-case[,-c(1:4)]

PAF$iska<-paste(PAF$ICD10,PAF$sex,PAF$kind,PAF$age)

data <- PAF %>%
  left_join(case, by = "iska")

data<-data[,-c(10)]

data_intinal <- data

# 找出包含NA的行
rows_with_na <- data[!complete.cases(data), ]
#print(rows_with_na)

#修正SH_smoking的case
#Smoking的P和PAF
exposure <- read_excel("~/Desktop/数据/广东省肿瘤归因/暴露/广东暴露率及RR（sex）.xlsx", sheet = "暴露率-sex+kind")
exposure <- exposure %>%
  filter(exposure == "Smoking")
exposure$ska<-paste(exposure$sex,exposure$kind,exposure$age)
exposure<-exposure[,-c(1:4)]

data1 <- data %>%
  filter(exposure == "Smoking" & ICD10 == "C33-C34")
data1$ska<-paste(data1$sex,data1$kind,data1$age)

data1 <- data1 %>%
  left_join(exposure, by = "ska")

colnames(data1)[colnames(data1) == "PAF"] <- "PAF_smoking"
data1<-data1[,c(7,12,13)]

#case_all(Second-hand smoking) = case_all(Second-hand smoking) *(1-PAF) * (1-P)
data2 <- data %>%
  filter(exposure == "Second-hand smoking" & ICD10 == "C33-C34")
data2$ska<-paste(data2$sex,data2$kind,data2$age)
data2 <- data2 %>%
  left_join(data1, by = "ska")

data2$inc_all<-data2$inc_all * (1-data2$rate) * (1-data2$PAF_smoking)
data2$mort_all<-data2$mort_all * (1-data2$rate) * (1-data2$PAF_smoking)

data2<-data2[,-c(12:14)]

data <- data %>%
  filter(exposure != "Second-hand smoking")

data <- rbind(data,data2)


#计算归因病例
data$PAF <- as.numeric(data$PAF)
data$PAF_lower <- as.numeric(data$PAF_lower)
data$PAF_upper <- as.numeric(data$PAF_upper)

data$inc_attribute<-data$inc_all*data$PAF
data$inc_attribute_lower<-data$inc_all*data$PAF_lower
data$inc_attribute_upper<-data$inc_all*data$PAF_upper
data$mort_attribute<-data$mort_all*data$PAF
data$mort_attribute_lower<-data$mort_all*data$PAF_lower
data$mort_attribute_upper<-data$mort_all*data$PAF_upper


#计算PAF（Second-hand smoking在全人群中的PAF）
data1 <- data_intinal %>%
  filter(exposure == "Smoking" & ICD10 == "C33-C34")
data1$ska<-paste(data1$sex,data1$kind,data1$age)
data1<-data1[,c(10:12)]

data2 <- data %>%
  filter(exposure == "Second-hand smoking" & ICD10 == "C33-C34")
data2$ska<-paste(data2$sex,data2$kind,data2$age)
data2<-data2[,-c(10,11)]
data2 <- data2 %>%
  left_join(data1, by = "ska")
data2<-data2[,-c(16)]
data2$PAF<-data2$inc_attribute/data2$inc_all
data2$PAF_lower<-data2$inc_attribute_lower/data2$inc_all
data2$PAF_upper<-data2$inc_attribute_upper/data2$inc_all

#names(data2)[names(data2) == "old_name"] <- "new_name"

data <- data %>%
  filter(exposure != "Second-hand smoking")

data <- rbind(data,data2)
data<-data[,-c(7:9)]

data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all


#合并rectum和colon
data$ICD10[data$ICD10 %in% c("C18","C19-C20")]<-"C18-C20"
data1 <- data %>%
  filter(ICD10 == "C18-C20")
data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+sex++kind+age+ICD10,data1,sum)
data1 <- data1 %>%
  mutate(site = "Colorectum")

data2 <- data1 %>%
  filter(exposure == "Alcohol")
data2$ska <- paste(data2$sex,data2$kind,data2$age)
data2<-data2[,c(6,7,15)]
data3 <- data1 %>%
  filter(exposure == "Physical inactivity")
data3$ska <- paste(data3$sex,data3$kind,data3$age)
data3<-data3[,-c(6,7)]
data3 <-data3%>%
  left_join(data2,by = "ska")
data3<-data3[,-c(13)]

data1 <- data1 %>%
  filter(exposure != "Physical inactivity")
data1<-rbind(data1,data3)

data1$PAF_inc<-data1$inc_attribute/data1$inc_all
data1$PAF_lower_inc<-data1$inc_attribute_lower/data1$inc_all
data1$PAF_upper_inc<-data1$inc_attribute_upper/data1$inc_all
data1$PAF_mort<-data1$mort_attribute/data1$mort_all
data1$PAF_lower_mort<-data1$mort_attribute_lower/data1$mort_all
data1$PAF_upper_mort<-data1$mort_attribute_upper/data1$mort_all

data <- data %>%
  filter(ICD10 != "C18-C20")
data<-rbind(data,data1)

#将Clonorchis sinensis的肝癌病例换成胆管癌
data1 <- data %>%
  filter(exposure == "Clonorchis sinensis" & ICD10 == "C22")
data1$inc_all1<-data1$inc_all * 0.125
data1$mort_all1<-data1$mort_all * 0.125

data1$inc_attribute<-data1$inc_all1*data1$PAF_inc
data1$inc_attribute_lower<-data1$inc_all1*data1$PAF_lower_inc
data1$inc_attribute_upper<-data1$inc_all1*data1$PAF_upper_inc
data1$mort_attribute<-data1$mort_all1*data1$PAF_mort
data1$mort_attribute_lower<-data1$mort_all1*data1$PAF_lower_mort
data1$mort_attribute_upper<-data1$mort_all1*data1$PAF_upper_mort

data1$PAF_inc<-data1$inc_attribute/data1$inc_all
data1$PAF_lower_inc<-data1$inc_attribute_lower/data1$inc_all
data1$PAF_upper_inc<-data1$inc_attribute_upper/data1$inc_all
data1$PAF_mort<-data1$mort_attribute/data1$mort_all
data1$PAF_lower_mort<-data1$mort_attribute_lower/data1$mort_all
data1$PAF_upper_mort<-data1$mort_attribute_upper/data1$mort_all
data1<-data1[,-c(21,22)]
data <- data %>%
  anti_join(data1, by = c("exposure", "ICD10"))
data<-rbind(data,data1)

#
data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1 <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, "(", data$inc_attribute_lower1, ",", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, "(", data$mort_attribute_lower1, ",", data$mort_attribute_upper1, ")", sep="")

data$PAF_inc2<-data$PAF_inc*100
data$PAF_lower_inc2<-data$PAF_lower_inc*100
data$PAF_upper_inc2<-data$PAF_upper_inc*100
data$PAF_mort2<-data$PAF_mort*100
data$PAF_lower_mort2<-data$PAF_lower_mort*100
data$PAF_upper_mort2<-data$PAF_upper_mort*100

data$PAF_inc2 <- sprintf("%0.2f", data$PAF_inc2)
data$PAF_lower_inc2 <- sprintf("%0.2f", data$PAF_lower_inc2)
data$PAF_upper_inc2 <- sprintf("%0.2f", data$PAF_upper_inc2)
data$PAF_95CI_inc <- paste(data$PAF_inc2, "(", data$PAF_lower_inc2, ",", data$PAF_upper_inc2, ")", sep="")

data$PAF_mort2 <- sprintf("%0.2f", data$PAF_mort2)
data$PAF_lower_mort2 <- sprintf("%0.2f", data$PAF_lower_mort2)
data$PAF_upper_mort2 <- sprintf("%0.2f", data$PAF_upper_mort2)
data$PAF_95CI_mort <- paste(data$PAF_mort2, "(", data$PAF_lower_mort2, ",", data$PAF_upper_mort2, ")", sep="")

data<-data[,-c(21:26, 29:34)]

na_rows <- data[is.na(data$PAF_inc), ]

library(openxlsx)
write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_age_case.xlsx")


#======================风险对的归因病例及PAF（分年龄、分城乡012、分性别012）============================================

data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_age_case.xlsx")

data<-data[,c(1:14)]

data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+sex+ICD10+site+age,data,sum)
data1 <- data1 %>%
  mutate(kind = 0)
data<-rbind(data,data1)

data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+kind+ICD10+site+age,data,sum)
data1 <- data1 %>%
  mutate(sex = 0)
data<-rbind(data,data1)

data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc2<-data$PAF_inc*100
data$PAF_lower_inc2<-data$PAF_lower_inc*100
data$PAF_upper_inc2<-data$PAF_upper_inc*100
data$PAF_mort2<-data$PAF_mort*100
data$PAF_lower_mort2<-data$PAF_lower_mort*100
data$PAF_upper_mort2<-data$PAF_upper_mort*100

data$PAF_inc2 <- sprintf("%0.2f", data$PAF_inc2)
data$PAF_lower_inc2 <- sprintf("%0.2f", data$PAF_lower_inc2)
data$PAF_upper_inc2 <- sprintf("%0.2f", data$PAF_upper_inc2)
data$PAF_95CI_inc <- paste(data$PAF_inc2, "(", data$PAF_lower_inc2, ",", data$PAF_upper_inc2, ")", sep="")

data$PAF_mort2 <- sprintf("%0.2f", data$PAF_mort2)
data$PAF_lower_mort2 <- sprintf("%0.2f", data$PAF_lower_mort2)
data$PAF_upper_mort2 <- sprintf("%0.2f", data$PAF_upper_mort2)
data$PAF_95CI_mort <- paste(data$PAF_mort2, "(", data$PAF_lower_mort2, ",", data$PAF_upper_mort2, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1 <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, "(", data$inc_attribute_lower1, ",", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, "(", data$mort_attribute_lower1, ",", data$mort_attribute_upper1, ")", sep="")

data<-data[,-c(21:26, 29:34)]

rows_with_na <- data[!complete.cases(data), ]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_age_kind012_sex012.xlsx")


#======================风险对的归因病例及PAF（全年龄段、分城乡、分性别）============================================

data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_age_kind012_sex012.xlsx")

data<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+sex+kind+ICD10+site,data,sum)

data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc2<-data$PAF_inc*100
data$PAF_lower_inc2<-data$PAF_lower_inc*100
data$PAF_upper_inc2<-data$PAF_upper_inc*100
data$PAF_mort2<-data$PAF_mort*100
data$PAF_lower_mort2<-data$PAF_lower_mort*100
data$PAF_upper_mort2<-data$PAF_upper_mort*100

data$PAF_inc2 <- sprintf("%0.2f", data$PAF_inc2)
data$PAF_lower_inc2 <- sprintf("%0.2f", data$PAF_lower_inc2)
data$PAF_upper_inc2 <- sprintf("%0.2f", data$PAF_upper_inc2)
data$PAF_95CI_inc <- paste(data$PAF_inc2, "(", data$PAF_lower_inc2, ",", data$PAF_upper_inc2, ")", sep="")

data$PAF_mort2 <- sprintf("%0.2f", data$PAF_mort2)
data$PAF_lower_mort2 <- sprintf("%0.2f", data$PAF_lower_mort2)
data$PAF_upper_mort2 <- sprintf("%0.2f", data$PAF_upper_mort2)
data$PAF_95CI_mort <- paste(data$PAF_mort2, "(", data$PAF_lower_mort2, ",", data$PAF_upper_mort2, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1 <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, "(", data$inc_attribute_lower1, ",", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, "(", data$mort_attribute_lower1, ",", data$mort_attribute_upper1, ")", sep="")

data<-data[,-c(20:25, 28:33)]

rows_with_na <- data[!complete.cases(data), ]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_all_age.xlsx")



#=================所有危险因素对某一癌症的PAF（分年龄、分城乡、分性别）==========================================================

library(readxl)
library(dplyr)
library(tidyr)
library(boot)
library(purrr)


data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_age_kind012_sex012.xlsx")

#现代生活方式导致的PAF
data <- data %>% filter(exposure %in% c("Diabetes", "Excess bodyweight", "Low fruit intake","Low vegetable intake","Physical inactivity","Red meat"))



#糖尿病和BMI的PAF
data$kasi <- paste(data$kind,data$age,data$sex,data$ICD10)


# 找出每个 age 中 BMI 和 diabetes 的 PAF 较小的那一行
to_remove <- data %>%
  filter(exposure %in% c("Excess bodyweight", "Diabetes")) %>%
  group_by(kasi) %>%
  slice_min(order_by = PAF_inc, n = 1, with_ties = FALSE) %>%
  ungroup()

# 从原数据中移除这些行
result_inc <- anti_join(data, to_remove, by = c("kasi", "exposure", "PAF_inc"))

data <- result_inc

#data <- data %>% mutate(exposure = if_else(exposure == "Excess bodyweight", "Diabetes", exposure))


# 定义计算联合PAF的函数
calculate_joint_PAF_inc <- function(data) {1 - prod(1 - data$PAF_inc)}

# 对每个分组应用Bootstrap模拟
set.seed(123)  # 设置随机种子以确保可重复性
results_inc <- data %>%
  group_by(sex, kind, ICD10, age) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_inc(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10, age) %>%
  summarise(
    PAF_inc = calculate_joint_PAF_inc(data[[1]]),
    PAF_lower_inc = quantile(bootstraps, 0.025, na.rm = TRUE),
    PAF_upper_inc = quantile(bootstraps, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# 找出每个 age 中 BMI 和 diabetes 的 PAF 较小的那一行
to_remove <- data %>%
  filter(exposure %in% c("Excess bodyweight", "Diabetes")) %>%
  group_by(kasi) %>%
  slice_min(order_by = PAF_mort, n = 1, with_ties = FALSE) %>%
  ungroup()

# 从原数据中移除这些行
result_mort <- anti_join(data, to_remove, by = c("kasi", "exposure", "PAF_mort"))

data <- result_mort

# 定义计算联合PAF的函数
calculate_joint_PAF_mort <- function(data) {1 - prod(1 - data$PAF_mort)}

# 对每个分组应用Bootstrap模拟
set.seed(123)  # 设置随机种子以确保可重复性
results_mort <- data %>%
  group_by(sex, kind, ICD10, age) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_mort(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10, age) %>%
  summarise(
    PAF_mort = calculate_joint_PAF_mort(data[[1]]),
    PAF_lower_mort = quantile(bootstraps, 0.025, na.rm = TRUE),
    PAF_upper_mort = quantile(bootstraps, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


results_inc$Iska<-paste(results_inc$ICD10,results_inc$sex,results_inc$kind,results_inc$age)
results_mort$Iska<-paste(results_mort$ICD10,results_mort$sex,results_mort$kind,results_mort$age)
results_mort<-results_mort[,-c(1:4)]

data1 <- results_inc %>% 
  left_join(results_mort, by = "Iska")

#+site,inc_all,mort_all
data$Iska<-paste(data$ICD10,data$sex,data$kind,data$age)
data<-data[c(6,7,8,25)]
data <- distinct(data)

merged_data <- data1 %>%
  left_join(data, by = "Iska")

data<-merged_data
data<-data[,-c(8)]

rows_with_na <- data[!complete.cases(data), ]
# 删除包含NA值的行
#df <- na.omit(data)
#data<-df


data$PAF_inc2<-data$PAF_inc*100
data$PAF_lower_inc2<-data$PAF_lower_inc*100
data$PAF_upper_inc2<-data$PAF_upper_inc*100
data$PAF_mort2<-data$PAF_mort*100
data$PAF_lower_mort2<-data$PAF_lower_mort*100
data$PAF_upper_mort2<-data$PAF_upper_mort*100

data$PAF_inc2 <- sprintf("%0.2f", data$PAF_inc2)
data$PAF_lower_inc2 <- sprintf("%0.2f", data$PAF_lower_inc2)
data$PAF_upper_inc2 <- sprintf("%0.2f", data$PAF_upper_inc2)
data$PAF_95CI_inc <- paste(data$PAF_inc2, "(", data$PAF_lower_inc2, ",", data$PAF_upper_inc2, ")", sep="")

data$PAF_mort2 <- sprintf("%0.2f", data$PAF_mort2)
data$PAF_lower_mort2 <- sprintf("%0.2f", data$PAF_lower_mort2)
data$PAF_upper_mort2 <- sprintf("%0.2f", data$PAF_upper_mort2)
data$PAF_95CI_mort <- paste(data$PAF_mort2, "(", data$PAF_lower_mort2, ",", data$PAF_upper_mort2, ")", sep="")

data<-data[,-c(14:19)]

rows_with_na <- data[!complete.cases(data), ]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_all_exposure_ICD10.xlsx")

#=================================所有危险因素对某一癌症的归因病例（分年龄、分城乡、分性别）=========================


data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_all_exposure_ICD10.xlsx")

# 找出包含NA的行
rows_with_na <- data[!complete.cases(data), ]
print(rows_with_na)

data$PAF_inc <- as.numeric(data$PAF_inc)
data$PAF_lower_inc <- as.numeric(data$PAF_lower_inc)
data$PAF_upper_inc <- as.numeric(data$PAF_upper_inc)

data$inc_attribute<-data$inc_all*data$PAF_inc
data$inc_attribute_lower<-data$inc_all*data$PAF_lower_inc
data$inc_attribute_upper<-data$inc_all*data$PAF_upper_inc

data$PAF_mort <- as.numeric(data$PAF_mort)
data$PAF_lower_mort <- as.numeric(data$PAF_lower_mort)
data$PAF_upper_mort <- as.numeric(data$PAF_upper_mort)

data$mort_attribute<-data$mort_all*data$PAF_mort
data$mort_attribute_lower<-data$mort_all*data$PAF_lower_mort
data$mort_attribute_upper<-data$mort_all*data$PAF_upper_mort

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1 <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, "(", data$inc_attribute_lower1, ",", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, "(", data$mort_attribute_lower1, ",", data$mort_attribute_upper1, ")", sep="")
data<-data[,-c(22:27)]


write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_case_all_exposure_ICD10.xlsx")


#=================所有危险因素对某一癌症的PAF（分城乡、分性别）==========================================================

data<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")

#现代生活方式导致的PAF
data <- data %>% filter(exposure %in% c("Diabetes", "Excess bodyweight", "Low fruit intake","Low vegetable intake","Physical inactivity","Red meat"))



#糖尿病和BMI的PAF
data$ksi <- paste(data$kind,data$sex,data$ICD10)


# 找出每个 age 中 BMI 和 diabetes 的 PAF 较小的那一行
to_remove <- data %>%
  filter(exposure %in% c("Excess bodyweight", "Diabetes")) %>%
  group_by(ksi) %>%
  slice_min(order_by = PAF_inc, n = 1, with_ties = FALSE) %>%
  ungroup()

# 从原数据中移除这些行
result_inc <- anti_join(data, to_remove, by = c("ksi", "exposure", "PAF_inc"))

data <- result_inc

#data <- data %>% mutate(exposure = if_else(exposure == "Excess bodyweight", "Diabetes", exposure))


# 定义计算联合PAF的函数
calculate_joint_PAF_inc <- function(data) {1 - prod(1 - data$PAF_inc)}


# 对每个分组应用Bootstrap模拟
set.seed(123)  # 设置随机种子以确保可重复性
results_inc <- data %>%
  group_by(sex, kind, ICD10) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_inc(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10) %>%
  summarise(
    PAF_inc = calculate_joint_PAF_inc(data[[1]]),
    PAF_lower_inc = quantile(bootstraps, 0.025),
    PAF_upper_inc = quantile(bootstraps, 0.975),
    .groups = "drop"
  )


data<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")

#糖尿病和BMI的PAF
data$ksi <- paste(data$kind,data$sex,data$ICD10)


# 找出每个 age 中 BMI 和 diabetes 的 PAF 较小的那一行
to_remove <- data %>%
  filter(exposure %in% c("Excess bodyweight", "Diabetes")) %>%
  group_by(ksi) %>%
  slice_min(order_by = PAF_mort, n = 1, with_ties = FALSE) %>%
  ungroup()

# 从原数据中移除这些行
result_mort <- anti_join(data, to_remove, by = c("ksi", "exposure", "PAF_mort"))

data <- result_mort

# 定义计算联合PAF的函数
calculate_joint_PAF_mort <- function(data) {1 - prod(1 - data$PAF_mort)}

# 对每个分组应用Bootstrap模拟
set.seed(123)  # 设置随机种子以确保可重复性
results_mort <- data %>%
  group_by(sex, kind, ICD10) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_mort(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10) %>%
  summarise(
    PAF_mort = calculate_joint_PAF_mort(data[[1]]),
    PAF_lower_mort = quantile(bootstraps, 0.025),
    PAF_upper_mort = quantile(bootstraps, 0.975),
    .groups = "drop"
  )


results_inc$Isk<-paste(results_inc$ICD10,results_inc$sex,results_inc$kind)
results_mort$Isk<-paste(results_mort$ICD10,results_mort$sex,results_mort$kind)
results_mort<-results_mort[,-c(1:3)]

data1 <- results_inc %>% 
  left_join(results_mort, by = "Isk")

#+site,inc_all,mort_all
data<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")
data$Isk<-paste(data$ICD10,data$sex,data$kind)
data<-data[c(5,6,7,24)]
data <- distinct(data)

merged_data <- data1 %>%
  left_join(data, by = "Isk")

data<-merged_data
data<-data[,-c(7)]


rows_with_na <- data[!complete.cases(data), ]
# 删除包含NA值的行
#df <- na.omit(data)
#data<-df

data$PAF_inc2<-data$PAF_inc*100
data$PAF_lower_inc2<-data$PAF_lower_inc*100
data$PAF_upper_inc2<-data$PAF_upper_inc*100
data$PAF_mort2<-data$PAF_mort*100
data$PAF_lower_mort2<-data$PAF_lower_mort*100
data$PAF_upper_mort2<-data$PAF_upper_mort*100

data$PAF_inc2 <- sprintf("%0.2f", data$PAF_inc2)
data$PAF_lower_inc2 <- sprintf("%0.2f", data$PAF_lower_inc2)
data$PAF_upper_inc2 <- sprintf("%0.2f", data$PAF_upper_inc2)
data$PAF_95CI_inc <- paste(data$PAF_inc2, "(", data$PAF_lower_inc2, ",", data$PAF_upper_inc2, ")", sep="")

data$PAF_mort2 <- sprintf("%0.2f", data$PAF_mort2)
data$PAF_lower_mort2 <- sprintf("%0.2f", data$PAF_lower_mort2)
data$PAF_upper_mort2 <- sprintf("%0.2f", data$PAF_upper_mort2)
data$PAF_95CI_mort <- paste(data$PAF_mort2, "(", data$PAF_lower_mort2, ",", data$PAF_upper_mort2, ")", sep="")

data<-data[,-c(13:18)]

data1 <- data[!is.na(data$PAF_inc), ]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_all_exposure_ICD10(all_age).xlsx")

#=================================所有危险因素对某一癌症的归因病例（分城乡、分性别）=========================
library(readxl)
library(dplyr)
library(tidyr)
library(boot)
library(purrr)
library(tidyverse)

data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_all_exposure_ICD10(all_age).xlsx")

data$PAF_inc <- as.numeric(data$PAF_inc)
data$PAF_lower_inc <- as.numeric(data$PAF_lower_inc)
data$PAF_upper_inc <- as.numeric(data$PAF_upper_inc)

data$inc_attribute<-data$inc_all*data$PAF_inc
data$inc_attribute_lower<-data$inc_all*data$PAF_lower_inc
data$inc_attribute_upper<-data$inc_all*data$PAF_upper_inc

data$PAF_mort <- as.numeric(data$PAF_mort)
data$PAF_lower_mort <- as.numeric(data$PAF_lower_mort)
data$PAF_upper_mort <- as.numeric(data$PAF_upper_mort)

data$mort_attribute<-data$mort_all*data$PAF_mort
data$mort_attribute_lower<-data$mort_all*data$PAF_lower_mort
data$mort_attribute_upper<-data$mort_all*data$PAF_upper_mort

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1 <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, "(", data$inc_attribute_lower1, ",", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, "(", data$mort_attribute_lower1, ",", data$mort_attribute_upper1, ")", sep="")
data<-data[,-c(21:26)]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_case_all_exposure_ICD10(all_age).xlsx")


#================================所有因素对所有癌症的归因病例和PAF（分年龄、分城乡、分性别）================================

case <- read_csv("~/Desktop/数据/广东省肿瘤归因/分析结果/cancer_all.csv")
case<-case[,-c(1,2,8)]

case <- case %>%
  filter(ICD10 == "ALL")

case <- case %>%
  filter(age != "f1" & age != "f2" & age != "f3" & age != "f4") 
case$age[case$age %in% c("f5","f6")]<-"20-29"
case$age[case$age %in% c("f7","f8")]<-"30-39"
case$age[case$age %in% c("f9","f10")]<-"40-49"
case$age[case$age %in% c("f11","f12")]<-"50-59"
case$age[case$age %in% c("f13","f14")]<-"60-69"
case$age[case$age %in% c("f15","f16","f17","f18")]<-"70-"
case<-aggregate(cbind(inc_all,mort_all)~sex+kind+ICD10,case,sum)
case$ska<-paste(case$sex,case$kind,case$age)
case<-case[,-c(1:3)]

data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/attribute_case_all_exposure_ICD10.xlsx")

data<-aggregate(cbind(inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~sex+kind+age,data,sum)
data$ska<-paste(data$sex,data$kind,data$age)
data<-data%>%
  left_join(case,by="ska")
data<-data[,-c(10)]

data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc2<-data$PAF_inc*100
data$PAF_lower_inc2<-data$PAF_lower_inc*100
data$PAF_upper_inc2<-data$PAF_upper_inc*100
data$PAF_mort2<-data$PAF_mort*100
data$PAF_lower_mort2<-data$PAF_lower_mort*100
data$PAF_upper_mort2<-data$PAF_upper_mort*100

data$PAF_inc2 <- sprintf("%0.2f", data$PAF_inc2)
data$PAF_lower_inc2 <- sprintf("%0.2f", data$PAF_lower_inc2)
data$PAF_upper_inc2 <- sprintf("%0.2f", data$PAF_upper_inc2)
data$PAF_95CI_inc <- paste(data$PAF_inc2, "(", data$PAF_lower_inc2, ",", data$PAF_upper_inc2, ")", sep="")

data$PAF_mort2 <- sprintf("%0.2f", data$PAF_mort2)
data$PAF_lower_mort2 <- sprintf("%0.2f", data$PAF_lower_mort2)
data$PAF_upper_mort2 <- sprintf("%0.2f", data$PAF_upper_mort2)
data$PAF_95CI_mort <- paste(data$PAF_mort2, "(", data$PAF_lower_mort2, ",", data$PAF_upper_mort2, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1 <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, "(", data$inc_attribute_lower1, ",", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, "(", data$mort_attribute_lower1, ",", data$mort_attribute_upper1, ")", sep="")
data<-data[,-c(18:23, 26:31)]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/所有危险因素-所有癌症PAF及case(age_group).xlsx")


#================================所有因素对所有癌症的归因病例和PAF（分城乡、分性别）================================

data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/所有危险因素-所有癌症PAF及case(age_group).xlsx")

data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~sex+kind,data,sum)
data<-data1
data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc2<-data$PAF_inc*100
data$PAF_lower_inc2<-data$PAF_lower_inc*100
data$PAF_upper_inc2<-data$PAF_upper_inc*100
data$PAF_mort2<-data$PAF_mort*100
data$PAF_lower_mort2<-data$PAF_lower_mort*100
data$PAF_upper_mort2<-data$PAF_upper_mort*100

data$PAF_inc2 <- sprintf("%0.2f", data$PAF_inc2)
data$PAF_lower_inc2 <- sprintf("%0.2f", data$PAF_lower_inc2)
data$PAF_upper_inc2 <- sprintf("%0.2f", data$PAF_upper_inc2)
data$PAF_95CI_inc <- paste(data$PAF_inc2, "(", data$PAF_lower_inc2, ",", data$PAF_upper_inc2, ")", sep="")

data$PAF_mort2 <- sprintf("%0.2f", data$PAF_mort2)
data$PAF_lower_mort2 <- sprintf("%0.2f", data$PAF_lower_mort2)
data$PAF_upper_mort2 <- sprintf("%0.2f", data$PAF_upper_mort2)
data$PAF_95CI_mort <- paste(data$PAF_mort2, "(", data$PAF_lower_mort2, ",", data$PAF_upper_mort2, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1 <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, "(", data$inc_attribute_lower1, ",", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, "(", data$mort_attribute_lower1, ",", data$mort_attribute_upper1, ")", sep="")
data<-data[,-c(17:22, 25:30)]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/所有危险因素-所有癌症PAF及case.xlsx")

#=================================某一危险因素对所有癌症的归因病例和PAF（分城乡、分性别）=========================

data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_all_age.xlsx")
data<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+sex+kind,data,sum)

data$sk<-paste(data$sex,data$kind)
data<-data[,-c(4,5)]

library(readr)
case <- read_csv("~/Desktop/数据/广东省肿瘤归因/cancer_all.csv")
case<-case[,-c(1,2,8)]

case <- case %>%
  filter(age != "f1" & age != "f2" & age != "f3" & age != "f4") 

#case <- case %>%
# filter(ICD10 == "C00-C14,except C11" | ICD10 == "C11" | ICD10 == "C15" | ICD10 == "C16" |ICD10 == "C18" | 
#         ICD10 == "C19-C20" | ICD10 == "C22" | ICD10 == "C23-C24" | ICD10 == "C25" | ICD10 == "C32" |
#        ICD10 == "C33-C34" | ICD10 == "C50" | ICD10 == "C51-C52" | ICD10 == "C53" | ICD10 == "C54" | ICD10 == "C56" |
#       ICD10 == "C60" | ICD10 == "C61" | ICD10 == "C64" | ICD10 == "C67" | ICD10 == "C73" | ICD10 == "C81" | ICD10 == "C92-C94")
case <- case %>%
  filter(ICD10 == "ALL")

case1<-aggregate(cbind(inc_all,mort_all)~sex+kind,case,sum)
case<-case1
case$sk<-paste(case$sex,case$kind)
case<-case[,-c(1,2)]

data <- data %>%
  left_join(case, by = "sk")
data<-data[,-c(10)]


data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc2<-data$PAF_inc*100
data$PAF_lower_inc2<-data$PAF_lower_inc*100
data$PAF_upper_inc2<-data$PAF_upper_inc*100
data$PAF_mort2<-data$PAF_mort*100
data$PAF_lower_mort2<-data$PAF_lower_mort*100
data$PAF_upper_mort2<-data$PAF_upper_mort*100

data$PAF_inc2 <- sprintf("%0.2f", data$PAF_inc2)
data$PAF_lower_inc2 <- sprintf("%0.2f", data$PAF_lower_inc2)
data$PAF_upper_inc2 <- sprintf("%0.2f", data$PAF_upper_inc2)
data$PAF_95CI_inc <- paste(data$PAF_inc2, "(", data$PAF_lower_inc2, ",", data$PAF_upper_inc2, ")", sep="")

data$PAF_mort2 <- sprintf("%0.2f", data$PAF_mort2)
data$PAF_lower_mort2 <- sprintf("%0.2f", data$PAF_lower_mort2)
data$PAF_upper_mort2 <- sprintf("%0.2f", data$PAF_upper_mort2)
data$PAF_95CI_mort <- paste(data$PAF_mort2, "(", data$PAF_lower_mort2, ",", data$PAF_upper_mort2, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1 <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, "(", data$inc_attribute_lower1, ",", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, "(", data$mort_attribute_lower1, ",", data$mort_attribute_upper1, ")", sep="")
data<-data[,-c(18:23, 26:31)]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/某一危险因素所有癌症PAF(所有癌种)及case.xlsx")

#=================================某一危险因素对所有癌症的归因病例和PAF（分年龄段、分城乡、分性别）=========================

data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_age_kind012_sex012.xlsx")
data<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+sex+kind+age,data,sum)

data$ska<-paste(data$sex,data$kind,data$age)
data<-data[,-c(5,6)]

library(readr)
case <- read_csv("~/Desktop/数据/广东省肿瘤归因/分析结果/cancer_all.csv")
case<-case[,-c(1,2,8)]

case <- case %>%
  filter(age != "f1" & age != "f2" & age != "f3" & age != "f4") 
case$age[case$age %in% c("f5","f6")]<-"20-29"
case$age[case$age %in% c("f7","f8")]<-"30-39"
case$age[case$age %in% c("f9","f10")]<-"40-49"
case$age[case$age %in% c("f11","f12")]<-"50-59"
case$age[case$age %in% c("f13","f14")]<-"60-69"
case$age[case$age %in% c("f15","f16","f17","f18")]<-"70-"

case <- case %>%
  filter(ICD10 == "ALL")

case1<-aggregate(cbind(inc_all,mort_all)~sex+kind+age,case,sum)
case<-case1
case$ska<-paste(case$sex,case$kind,case$age)
case<-case[,-c(1:3)]

data <- data %>%
  left_join(case, by = "ska")
data<-data[,-c(11)]


data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc2<-data$PAF_inc*100
data$PAF_lower_inc2<-data$PAF_lower_inc*100
data$PAF_upper_inc2<-data$PAF_upper_inc*100
data$PAF_mort2<-data$PAF_mort*100
data$PAF_lower_mort2<-data$PAF_lower_mort*100
data$PAF_upper_mort2<-data$PAF_upper_mort*100

data$PAF_inc2 <- sprintf("%0.2f", data$PAF_inc2)
data$PAF_lower_inc2 <- sprintf("%0.2f", data$PAF_lower_inc2)
data$PAF_upper_inc2 <- sprintf("%0.2f", data$PAF_upper_inc2)
data$PAF_95CI_inc <- paste(data$PAF_inc2, "(", data$PAF_lower_inc2, ",", data$PAF_upper_inc2, ")", sep="")

data$PAF_mort2 <- sprintf("%0.2f", data$PAF_mort2)
data$PAF_lower_mort2 <- sprintf("%0.2f", data$PAF_lower_mort2)
data$PAF_upper_mort2 <- sprintf("%0.2f", data$PAF_upper_mort2)
data$PAF_95CI_mort <- paste(data$PAF_mort2, "(", data$PAF_lower_mort2, ",", data$PAF_upper_mort2, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1 <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, "(", data$inc_attribute_lower1, ",", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, "(", data$mort_attribute_lower1, ",", data$mort_attribute_upper1, ")", sep="")
data<-data[,-c(19:24, 27:32)]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/某一危险因素所有癌症PAF(所有癌种)及case_divided_age.xlsx")


#=================================某一factor内所有危险因素对所有癌症的归因病例（分城乡、分性别）=========================
library(readxl)
library(dplyr)
library(tidyr)

#连上factor
data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")
data1 <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/factor.xlsx")

#data1 <- data1[,c(1,2)]
#data1 <- distinct(data1)
data <- data %>%
  left_join(data1, by = "exposure")


#data1 <- data %>%
#filter(exposure == "Physical inactivity"|exposure == "Diabetes"|exposure == "Excess bodyweight"|exposure == "Red meat"|exposure == "Low fruit intake"|exposure == "Low vegetable intake")

data1 <- data %>%
  filter(factor == "Behaviour factor")
data2 <- data %>%
  filter(factor == "Metabolic factor")
data3 <- data %>%
  filter(factor == "Dietary factor")
data4 <- data %>%
  filter(factor == "Infectious agent")


#============分tobacco,dietary,infection================

data1 <- data %>% filter(factor == "Infectious factor")
data2 <- data %>% filter(factor == "Dietary factor")
data3 <- data %>% filter(factor == "Tobacco consumption")

#============分三个factor================
data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")

data$factor <- NA  # 创建一个空的factor列

data$factor[data$exposure %in% c("Alcohol", "Smoking", "Second-hand smoking")] <- "Unhealthy behavior-related risk factors"
data$factor[data$exposure %in% c("Excess bodyweight", "Low vegetable intake", "Physical inactivity", "Diabetes", "Red meat", "Low fruit intake")] <- "Social development-related risk factors"
data$factor[data$exposure %in% c("Clonorchis sinensis", "Hepatitis C virus", "Helicobacter pylori", "Epstein-Barr virus", "Human papillomavirus", "Hepatitis B virus")] <- "Infection-related risk factors"

data1 <- data %>%
  filter(factor == "Unhealthy behavior-related risk factors")
data2 <- data %>%
  filter(factor == "Social development-related risk factors")
data3 <- data %>%
  filter(factor == "Infection-related risk factors")

#============================
#每个factor的单个癌种-所有危险因素

#Behaviour factor
calculate_joint_PAF_inc <- function(data) {1 - prod(1 - data$PAF_inc)}

set.seed(123) 
results_inc <- data1 %>%
  group_by(sex, kind, ICD10) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_inc(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10) %>%
  summarise(
    PAF_inc = calculate_joint_PAF_inc(data[[1]]),
    PAF_lower_inc = quantile(bootstraps, 0.025),
    PAF_upper_inc = quantile(bootstraps, 0.975),
    .groups = "drop"
  )

calculate_joint_PAF_mort <- function(data) {1 - prod(1 - data$PAF_mort)}

set.seed(123)
results_mort <- data1 %>%
  group_by(sex, kind, ICD10) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_mort(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10) %>%
  summarise(
    PAF_mort = calculate_joint_PAF_mort(data[[1]]),
    PAF_lower_mort = quantile(bootstraps, 0.025),
    PAF_upper_mort = quantile(bootstraps, 0.975),
    .groups = "drop"
  )


results_inc$Isk<-paste(results_inc$ICD10,results_inc$sex,results_inc$kind)
results_mort$Isk<-paste(results_mort$ICD10,results_mort$sex,results_mort$kind)
results_mort<-results_mort[,-c(1:3)]

data11 <- results_inc %>% 
  left_join(results_mort, by = "Isk")

#+site,inc_all,mort_all
data1$Isk<-paste(data1$ICD10,data1$sex,data1$kind)
data1<-data1[c(5,6,7,24,25)]
data1 <- distinct(data1)

merged_data1 <- data11 %>%
  left_join(data1, by = "Isk")

data1<-merged_data1
data1<-data1[,-c(7)]

#Metabolic factor
calculate_joint_PAF_inc <- function(data) {1 - prod(1 - data$PAF_inc)}

set.seed(123) 
results_inc <- data2 %>%
  group_by(sex, kind, ICD10) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_inc(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10) %>%
  summarise(
    PAF_inc = calculate_joint_PAF_inc(data[[1]]),
    PAF_lower_inc = quantile(bootstraps, 0.025),
    PAF_upper_inc = quantile(bootstraps, 0.975),
    .groups = "drop"
  )

calculate_joint_PAF_mort <- function(data) {1 - prod(1 - data$PAF_mort)}

set.seed(123)
results_mort <- data2 %>%
  group_by(sex, kind, ICD10) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_mort(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10) %>%
  summarise(
    PAF_mort = calculate_joint_PAF_mort(data[[1]]),
    PAF_lower_mort = quantile(bootstraps, 0.025),
    PAF_upper_mort = quantile(bootstraps, 0.975),
    .groups = "drop"
  )


results_inc$Isk<-paste(results_inc$ICD10,results_inc$sex,results_inc$kind)
results_mort$Isk<-paste(results_mort$ICD10,results_mort$sex,results_mort$kind)
results_mort<-results_mort[,-c(1:3)]

data21 <- results_inc %>% 
  left_join(results_mort, by = "Isk")

#+site,inc_all,mort_all
data2$Isk<-paste(data2$ICD10,data2$sex,data2$kind)
data2<-data2[c(5,6,7,24,25)]
data2 <- distinct(data2)

merged_data2 <- data21 %>%
  left_join(data2, by = "Isk")

data2<-merged_data2
data2<-data2[,-c(7)]


#Dietary factor
calculate_joint_PAF_inc <- function(data) {1 - prod(1 - data$PAF_inc)}

set.seed(123) 
results_inc <- data3 %>%
  group_by(sex, kind, ICD10) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_inc(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10) %>%
  summarise(
    PAF_inc = calculate_joint_PAF_inc(data[[1]]),
    PAF_lower_inc = quantile(bootstraps, 0.025),
    PAF_upper_inc = quantile(bootstraps, 0.975),
    .groups = "drop"
  )

calculate_joint_PAF_mort <- function(data) {1 - prod(1 - data$PAF_mort)}

set.seed(123)
results_mort <- data3 %>%
  group_by(sex, kind, ICD10) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_mort(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10) %>%
  summarise(
    PAF_mort = calculate_joint_PAF_mort(data[[1]]),
    PAF_lower_mort = quantile(bootstraps, 0.025),
    PAF_upper_mort = quantile(bootstraps, 0.975),
    .groups = "drop"
  )


results_inc$Isk<-paste(results_inc$ICD10,results_inc$sex,results_inc$kind)
results_mort$Isk<-paste(results_mort$ICD10,results_mort$sex,results_mort$kind)
results_mort<-results_mort[,-c(1:3)]

data31 <- results_inc %>% 
  left_join(results_mort, by = "Isk")

#+site,inc_all,mort_all
data3$Isk<-paste(data3$ICD10,data3$sex,data3$kind)
data3<-data3[c(5,6,7,24,25)]
data3 <- distinct(data3)

merged_data3 <- data31 %>%
  left_join(data3, by = "Isk")

data3<-merged_data3
data3<-data3[,-c(7)]

#Infectious agent
calculate_joint_PAF_inc <- function(data) {1 - prod(1 - data$PAF_inc)}

set.seed(123) 
results_inc <- data4 %>%
  group_by(sex, kind, ICD10) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_inc(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10) %>%
  summarise(
    PAF_inc = calculate_joint_PAF_inc(data[[1]]),
    PAF_lower_inc = quantile(bootstraps, 0.025),
    PAF_upper_inc = quantile(bootstraps, 0.975),
    .groups = "drop"
  )

calculate_joint_PAF_mort <- function(data) {1 - prod(1 - data$PAF_mort)}

set.seed(123)
results_mort <- data4 %>%
  group_by(sex, kind, ICD10) %>%
  nest() %>%
  mutate(bootstraps = list(replicate(5000, {
    sampled_data <- data[[1]][sample(nrow(data[[1]]), replace = TRUE), ]
    calculate_joint_PAF_mort(sampled_data)
  }))) %>%
  unnest(bootstraps) %>%
  group_by(sex, kind, ICD10) %>%
  summarise(
    PAF_mort = calculate_joint_PAF_mort(data[[1]]),
    PAF_lower_mort = quantile(bootstraps, 0.025),
    PAF_upper_mort = quantile(bootstraps, 0.975),
    .groups = "drop"
  )


results_inc$Isk<-paste(results_inc$ICD10,results_inc$sex,results_inc$kind)
results_mort$Isk<-paste(results_mort$ICD10,results_mort$sex,results_mort$kind)
results_mort<-results_mort[,-c(1:3)]

data41 <- results_inc %>% 
  left_join(results_mort, by = "Isk")

#+site,inc_all,mort_all
data4$Isk<-paste(data4$ICD10,data4$sex,data4$kind)
data4<-data4[c(5,6,7,24,25)]
data4 <- distinct(data4)

merged_data4 <- data41 %>%
  left_join(data4, by = "Isk")

data4<-merged_data4
data4<-data4[,-c(7)]

data<-rbind(data1,data2,data3,data4)


rows_with_na <- data[!complete.cases(data), ]
# 删除包含NA值的行
#df <- na.omit(data)
#data<-df

data$PAF_inc <- as.numeric(data$PAF_inc)
data$PAF_lower_inc <- as.numeric(data$PAF_lower_inc)
data$PAF_upper_inc <- as.numeric(data$PAF_upper_inc)

data$inc_attribute<-data$inc_all*data$PAF_inc
data$inc_attribute_lower<-data$inc_all*data$PAF_lower_inc
data$inc_attribute_upper<-data$inc_all*data$PAF_upper_inc

data$PAF_mort <- as.numeric(data$PAF_mort)
data$PAF_lower_mort <- as.numeric(data$PAF_lower_mort)
data$PAF_upper_mort <- as.numeric(data$PAF_upper_mort)

data$mort_attribute<-data$mort_all*data$PAF_mort
data$mort_attribute_lower<-data$mort_all*data$PAF_lower_mort
data$mort_attribute_upper<-data$mort_all*data$PAF_upper_mort

#合并所有癌种
data1<-aggregate(cbind(inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~sex+kind+factor,data,sum)
data<-data1
data$sk<-paste(data$sex,data$kind)

#连上全癌种的inc_all,mort_all
library(readr)
case <- read_csv("~/Desktop/数据/广东省肿瘤归因/分析结果/cancer_all.csv")
case<-case[,-c(1,2,8)]

case <- case %>%
  filter(age != "f1" & age != "f2" & age != "f3" & age != "f4") 

#case <- case %>%
# filter(ICD10 == "C00-C14,except C11" | ICD10 == "C11" | ICD10 == "C15" | ICD10 == "C16" |ICD10 == "C18" | 
#         ICD10 == "C19-C20" | ICD10 == "C22" | ICD10 == "C23-C24" | ICD10 == "C25" | ICD10 == "C32" |
#        ICD10 == "C33-C34" | ICD10 == "C50" | ICD10 == "C51-C52" | ICD10 == "C53" | ICD10 == "C54" | ICD10 == "C56" |
#       ICD10 == "C60" | ICD10 == "C61" | ICD10 == "C64" | ICD10 == "C67" | ICD10 == "C73" | ICD10 == "C81" | ICD10 == "C92-C94")

case <- case %>%
  filter(ICD10 == "ALL")

case1<-aggregate(cbind(inc_all,mort_all)~sex+kind,case,sum)
case<-case1
case$sk<-paste(case$sex,case$kind)
case<-case[,-c(1,2)]

data <- data %>%
  left_join(case, by = "sk")


data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc2<-data$PAF_inc*100
data$PAF_lower_inc2<-data$PAF_lower_inc*100
data$PAF_upper_inc2<-data$PAF_upper_inc*100
data$PAF_mort2<-data$PAF_mort*100
data$PAF_lower_mort2<-data$PAF_lower_mort*100
data$PAF_upper_mort2<-data$PAF_upper_mort*100

data$PAF_inc2 <- sprintf("%0.2f", data$PAF_inc2)
data$PAF_lower_inc2 <- sprintf("%0.2f", data$PAF_lower_inc2)
data$PAF_upper_inc2 <- sprintf("%0.2f", data$PAF_upper_inc2)
data$PAF_95CI_inc <- paste(data$PAF_inc2, "(", data$PAF_lower_inc2, ",", data$PAF_upper_inc2, ")", sep="")

data$PAF_mort2 <- sprintf("%0.2f", data$PAF_mort2)
data$PAF_lower_mort2 <- sprintf("%0.2f", data$PAF_lower_mort2)
data$PAF_upper_mort2 <- sprintf("%0.2f", data$PAF_upper_mort2)
data$PAF_95CI_mort <- paste(data$PAF_mort2, "(", data$PAF_lower_mort2, ",", data$PAF_upper_mort2, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1 <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, "(", data$inc_attribute_lower1, ",", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, "(", data$mort_attribute_lower1, ",", data$mort_attribute_upper1, ")", sep="")
data<-data[,-c(10, 19:24, 27:32)]

library(openxlsx)
write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/某一factor所有癌症PAF及case.xlsx")


#======================各个factor分癌种的PAF（在所有癌种中)===================================
#连上上面的data1，data2，data3，data4
data<-rbind(data1,data2,data3,data4)

data$PAF_inc <- as.numeric(data$PAF_inc)
data$PAF_lower_inc <- as.numeric(data$PAF_lower_inc)
data$PAF_upper_inc <- as.numeric(data$PAF_upper_inc)
data$inc_attribute<-data$inc_all*data$PAF_inc
data$inc_attribute_lower<-data$inc_all*data$PAF_lower_inc
data$inc_attribute_upper<-data$inc_all*data$PAF_upper_inc

data$PAF_mort <- as.numeric(data$PAF_mort)
data$PAF_lower_mort <- as.numeric(data$PAF_lower_mort)
data$PAF_upper_mort <- as.numeric(data$PAF_upper_mort)
data$mort_attribute<-data$mort_all*data$PAF_mort
data$mort_attribute_lower<-data$mort_all*data$PAF_lower_mort
data$mort_attribute_upper<-data$mort_all*data$PAF_upper_mort

data<-data[,-c(11,12)]
data$sk<-paste(data$sex,data$kind)

data <- data %>%
  left_join(case, by = "sk")
data<-data[,-c(18)]


data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc1<-data$PAF_inc*100
data$PAF_lower_inc1<-data$PAF_lower_inc*100
data$PAF_upper_inc1<-data$PAF_upper_inc*100
data$PAF_mort1<-data$PAF_mort*100
data$PAF_lower_mort1<-data$PAF_lower_mort*100
data$PAF_upper_mort1<-data$PAF_upper_mort*100

data$PAF_inc1 <- sprintf("%0.2f", data$PAF_inc1)
data$PAF_lower_inc1 <- sprintf("%0.2f", data$PAF_lower_inc1)
data$PAF_upper_inc1 <- sprintf("%0.2f", data$PAF_upper_inc1)
data$PAF_95CI_inc <- paste(data$PAF_inc1, " (", data$PAF_lower_inc1, ", ", data$PAF_upper_inc1, ")", sep="")

data$PAF_mort1 <- sprintf("%0.2f", data$PAF_mort1)
data$PAF_lower_mort1 <- sprintf("%0.2f", data$PAF_lower_mort1)
data$PAF_upper_mort1 <- sprintf("%0.2f", data$PAF_upper_mort1)
data$PAF_95CI_mort <- paste(data$PAF_mort1, " (", data$PAF_lower_mort1, ", ", data$PAF_upper_mort1, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1<- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, " (", data$inc_attribute_lower1, ", ", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, " (", data$mort_attribute_lower1, ", ", data$mort_attribute_upper1, ")", sep="")

data<-data[,-c(20:25, 28:33)]
write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/各个factor分癌种的PAF.xlsx")




data1<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/各个factor分癌种的PAF.xlsx")
data2<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/某一factor所有癌症PAF及case.xlsx")

data2 <- data2 %>%
  mutate(ICD10 = "All")
data2 <- data2 %>%
  mutate(site = "All")

data<-rbind(data1,data2)

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/各个factor分癌种+所有癌症的PAF.xlsx")


#=====================================删除男女合计和城乡合计，重新合起来===========================================
#1
data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_case_all_exposure_ICD10-unhealthy style.xlsx")

data <- data %>%
  filter(sex != 0)

data <- data %>%
  filter(kind != 0)

# 将所有 NA 值替换为 0
data <- data %>%
  mutate(across(everything(), ~ifelse(is.na(.), 0, .)))

data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~sex+age+ICD10+site,data,sum)
data1 <- data1 %>%
  mutate(kind = 0)
data<-data[,-c(5:10,14:15,22:23)]
data<-rbind(data,data1)

data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~kind+age+ICD10+site,data,sum)
data1 <- data1 %>%
  mutate(sex = 0)
data<-rbind(data,data1)
table(data$sex)
table(data$kind)


data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc1<-data$PAF_inc*100
data$PAF_lower_inc1<-data$PAF_lower_inc*100
data$PAF_upper_inc1<-data$PAF_upper_inc*100
data$PAF_mort1<-data$PAF_mort*100
data$PAF_lower_mort1<-data$PAF_lower_mort*100
data$PAF_upper_mort1<-data$PAF_upper_mort*100

data$PAF_inc1 <- sprintf("%0.2f", data$PAF_inc1)
data$PAF_lower_inc1 <- sprintf("%0.2f", data$PAF_lower_inc1)
data$PAF_upper_inc1 <- sprintf("%0.2f", data$PAF_upper_inc1)
data$PAF_95CI_inc <- paste(data$PAF_inc1, " (", data$PAF_lower_inc1, ", ", data$PAF_upper_inc1, ")", sep="")

data$PAF_mort1 <- sprintf("%0.2f", data$PAF_mort1)
data$PAF_lower_mort1 <- sprintf("%0.2f", data$PAF_lower_mort1)
data$PAF_upper_mort1 <- sprintf("%0.2f", data$PAF_upper_mort1)
data$PAF_95CI_mort <- paste(data$PAF_mort1, " (", data$PAF_lower_mort1, ", ", data$PAF_upper_mort1, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1<- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, " (", data$inc_attribute_lower1, ", ", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, " (", data$mort_attribute_lower1, ", ", data$mort_attribute_upper1, ")", sep="")

data<-data[,-c(20:25, 28:33)]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/attribute_case_all_exposure_ICD10-unhealthy style.xlsx")



#2
data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_case_all_exposure_ICD10(all_age).xlsx")

data <- data %>%
  filter(sex != 0)

data <- data %>%
  filter(kind != 0)

# 将所有 NA 值替换为 0
data <- data %>%
  mutate(across(everything(), ~ifelse(is.na(.), 0, .)))

data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~sex+ICD10+site,data,sum)
data1 <- data1 %>%
  mutate(kind = 0)
data<-data[,-c(4:9,13:14,21:22)]
data<-rbind(data,data1)

data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~kind+ICD10+site,data,sum)
data1 <- data1 %>%
  mutate(sex = 0)
data<-rbind(data,data1)
table(data$sex)
table(data$kind)


data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc1<-data$PAF_inc*100
data$PAF_lower_inc1<-data$PAF_lower_inc*100
data$PAF_upper_inc1<-data$PAF_upper_inc*100
data$PAF_mort1<-data$PAF_mort*100
data$PAF_lower_mort1<-data$PAF_lower_mort*100
data$PAF_upper_mort1<-data$PAF_upper_mort*100

data$PAF_inc1 <- sprintf("%0.2f", data$PAF_inc1)
data$PAF_lower_inc1 <- sprintf("%0.2f", data$PAF_lower_inc1)
data$PAF_upper_inc1 <- sprintf("%0.2f", data$PAF_upper_inc1)
data$PAF_95CI_inc <- paste(data$PAF_inc1, " (", data$PAF_lower_inc1, ", ", data$PAF_upper_inc1, ")", sep="")

data$PAF_mort1 <- sprintf("%0.2f", data$PAF_mort1)
data$PAF_lower_mort1 <- sprintf("%0.2f", data$PAF_lower_mort1)
data$PAF_upper_mort1 <- sprintf("%0.2f", data$PAF_upper_mort1)
data$PAF_95CI_mort <- paste(data$PAF_mort1, " (", data$PAF_lower_mort1, ", ", data$PAF_upper_mort1, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1<- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, " (", data$inc_attribute_lower1, ", ", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, " (", data$mort_attribute_lower1, ", ", data$mort_attribute_upper1, ")", sep="")

data<-data[,-c(19:24, 27:32)]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/attribute_case_all_exposure_ICD10(all_age).xlsx")

#3
data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/某一factor所有癌症PAF及case.xlsx")
data <- data %>%
  filter(sex != 0)

data <- data %>%
  filter(kind != 0)

# 将所有 NA 值替换为 0
data <- data %>%
  mutate(across(everything(), ~ifelse(is.na(.), 0, .)))

data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~sex+factor,data,sum)
data1 <- data1 %>%
  mutate(kind = 0)
data<-data[,-c(12:21)]
data<-rbind(data,data1)

data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~kind+factor,data,sum)
data1 <- data1 %>%
  mutate(sex = 0)
data<-rbind(data,data1)
table(data$sex)
table(data$kind)


data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc1<-data$PAF_inc*100
data$PAF_lower_inc1<-data$PAF_lower_inc*100
data$PAF_upper_inc1<-data$PAF_upper_inc*100
data$PAF_mort1<-data$PAF_mort*100
data$PAF_lower_mort1<-data$PAF_lower_mort*100
data$PAF_upper_mort1<-data$PAF_upper_mort*100

data$PAF_inc1 <- sprintf("%0.2f", data$PAF_inc1)
data$PAF_lower_inc1 <- sprintf("%0.2f", data$PAF_lower_inc1)
data$PAF_upper_inc1 <- sprintf("%0.2f", data$PAF_upper_inc1)
data$PAF_95CI_inc <- paste(data$PAF_inc1, " (", data$PAF_lower_inc1, ", ", data$PAF_upper_inc1, ")", sep="")

data$PAF_mort1 <- sprintf("%0.2f", data$PAF_mort1)
data$PAF_lower_mort1 <- sprintf("%0.2f", data$PAF_lower_mort1)
data$PAF_upper_mort1 <- sprintf("%0.2f", data$PAF_upper_mort1)
data$PAF_95CI_mort <- paste(data$PAF_mort1, " (", data$PAF_lower_mort1, ", ", data$PAF_upper_mort1, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1<- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, " (", data$inc_attribute_lower1, ", ", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, " (", data$mort_attribute_lower1, ", ", data$mort_attribute_upper1, ")", sep="")

data<-data[,-c(18:23, 26:31)]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/某一factor所有癌症PAF及case.xlsx")

#4

data <- read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/各个factor分癌种的PAF.xlsx")
data <- data %>%
  filter(sex != 0)

data <- data %>%
  filter(kind != 0)

# 将所有 NA 值替换为 0
data <- data %>%
  mutate(across(everything(), ~ifelse(is.na(.), 0, .)))

data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~sex+factor+ICD10+site,data,sum)
data1 <- data1 %>%
  mutate(kind = 0)
data<-data[,-c(4:9,20:23)]
data<-rbind(data,data1)

data1<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~kind+factor+ICD10+site,data,sum)
data1 <- data1 %>%
  mutate(sex = 0)
data<-rbind(data,data1)
table(data$sex)
table(data$kind)


data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all

data$PAF_inc1<-data$PAF_inc*100
data$PAF_lower_inc1<-data$PAF_lower_inc*100
data$PAF_upper_inc1<-data$PAF_upper_inc*100
data$PAF_mort1<-data$PAF_mort*100
data$PAF_lower_mort1<-data$PAF_lower_mort*100
data$PAF_upper_mort1<-data$PAF_upper_mort*100

data$PAF_inc1 <- sprintf("%0.2f", data$PAF_inc1)
data$PAF_lower_inc1 <- sprintf("%0.2f", data$PAF_lower_inc1)
data$PAF_upper_inc1 <- sprintf("%0.2f", data$PAF_upper_inc1)
data$PAF_95CI_inc <- paste(data$PAF_inc1, " (", data$PAF_lower_inc1, ", ", data$PAF_upper_inc1, ")", sep="")

data$PAF_mort1 <- sprintf("%0.2f", data$PAF_mort1)
data$PAF_lower_mort1 <- sprintf("%0.2f", data$PAF_lower_mort1)
data$PAF_upper_mort1 <- sprintf("%0.2f", data$PAF_upper_mort1)
data$PAF_95CI_mort <- paste(data$PAF_mort1, " (", data$PAF_lower_mort1, ", ", data$PAF_upper_mort1, ")", sep="")

data$mort_attribute1 <- as.integer(data$mort_attribute)
data$mort_attribute_lower1<- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper1 <- as.integer(data$mort_attribute_upper)
data$inc_attribute1 <- as.integer(data$inc_attribute)
data$inc_attribute_lower1 <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper1 <- as.integer(data$inc_attribute_upper)

data$inc_95CI <- paste(data$inc_attribute1, " (", data$inc_attribute_lower1, ", ", data$inc_attribute_upper1, ")", sep="")
data$mort_95CI <- paste(data$mort_attribute1, " (", data$mort_attribute_lower1, ", ", data$mort_attribute_upper1, ")", sep="")

data<-data[,-c(20:25, 28:33)]

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/各个factor分癌种的PAF.xlsx")



data1<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/各个factor分癌种的PAF.xlsx")
data2<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/某一factor所有癌症PAF及case.xlsx")

data2 <- data2 %>%
  mutate(ICD10 = "All")
data2 <- data2 %>%
  mutate(site = "All")

data<-rbind(data1,data2)

write.xlsx(data, "~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/各个factor分癌种+所有癌症的PAF.xlsx")


#======多分类exposure的PAF计算=================================

library(dplyr)
library(tidyr)

# 定义计算单个危险因素PAF的函数
calculate_single_PAF <- function(data) {
  numerator <- sum(data$P * (data$RR - 1))
  denominator <- numerator + 1
  PAF <- numerator / denominator
  return(PAF)
}

# 计算每个危险因素的PAF
paf_results <- data %>%
  group_by(kind) %>%
  summarise(PAF = calculate_single_PAF(cur_data_all()), .groups = "drop")

# Bootstrap模拟计算95%CI
set.seed(123) 
results_with_ci <- data %>%
  group_by(kind) %>%
  nest() %>%
  mutate(
    bootstraps = map(data, ~ replicate(5000, {
      sampled_data <- .x[sample(nrow(.x), replace = TRUE), ]
      calculate_single_PAF(sampled_data)
    }, simplify = TRUE))
  ) %>%
  unnest_wider(bootstraps) %>%
  mutate(
    PAF = map_dbl(data, calculate_single_PAF),
    PAF_lower = map_dbl(bootstraps, ~ quantile(.x, 0.025)),
    PAF_upper = map_dbl(bootstraps, ~ quantile(.x, 0.975))
  ) %>%
  select(kind, PAF, PAF_lower, PAF_upper)

print(results_with_ci)





#==================================最新：按癌种区分的所有危险因素PAF条图(仅有男女合计）====================================
library(gridExtra)
library(ggplot2)
library(grid)
library(readxl)
library(dplyr)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_case_all_exposure_ICD10(all_age).xlsx")

data1 <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/所有危险因素-所有癌症PAF及case.xlsx")
data1 <- data1 %>%
  mutate(site = "All")
data1 <- data1 %>%
  mutate(ICD10 = "All")

data<-rbind(data,data1)

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 0)

data$PAF_inc<-data$PAF_inc*100
data$PAF_lower_inc<-data$PAF_lower_inc*100
data$PAF_upper_inc<-data$PAF_upper_inc*100
data$PAF_inc1 <- sprintf("%0.1f", data$PAF_inc)

data$PAF_mort<-data$PAF_mort*100
data$PAF_lower_mort<-data$PAF_lower_mort*100
data$PAF_upper_mort<-data$PAF_upper_mort*100
data$PAF_mort1 <- sprintf("%0.1f", data$PAF_mort)

data1 <- data %>% filter(sex == 0)
data1 <- data1 %>%
  arrange(PAF_inc) %>%  # 对每个 sex 组内按 PAF_mort 排序
  mutate(site = factor(site, levels = unique(site))) # 确保排序后不改变 site 的因子级别

ggplot(data1, aes(x = PAF_inc, y = site)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, fill = "#3366CC") +  # 设置条形颜色
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100), labels = scales::percent_format(scale = 1), expand = c(0, 0)) +  # 扩展x轴范围
  scale_y_discrete() +  # 使用 discrete scale for y axis
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.text.y = element_text(size = 9, family = "Arial", color = "black"),  # 纵坐标标签样式
    axis.line.x = element_line(color = "grey", size = 0.3),
    axis.line.y = element_line(color = "grey", size = 0.3),
    axis.ticks = element_line(color = "grey", size = 0.3),
    legend.position = "bottom",  # 将图例放置在底部
    legend.title = element_blank(),  # 删除图例标题
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_text(size = 9, family = "Arial", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
    strip.placement = "left",  # 将facet标签放置在左侧
    plot.margin = margin(l = 50),  # 调整左边距以腾出空间给文本标签
    axis.text.x = element_text(family = "Arial", color = "black")  # 设置横坐标轴字体
  ) +
  # 添加inc_attribute文本标签到x = 115的位置
  geom_text(aes(x = 115 , y = site, label = inc_95CI), 
            hjust = 0, size = 2.5, family = "Arial", color = "black") +  # 将标签放到x = 115位置
  coord_cartesian(clip = "off") +
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = PAF_inc + 1, y = site, label = PAF_inc1), 
            vjust = 0.5, hjust = 0, size = 2.5, family = "Arial", color = "black") +
  # 添加"Attributable Cases"标签
  geom_text(data = data1, 
            aes(x = 150, y = 24.5, label = "Attributable Cases"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Arial", color = "black") + 
  # 添加"PAF(%)"标签
  geom_text(data = data1, 
            aes(x = 60, y = 24.5, label = "PAF(%)"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Arial", color = "black") +
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100),labels = c("0", "50", "100"), expand = c(0, 0))  # 扩展x轴范围并设置刻度标签



data2 <- data %>% filter(sex == 0)
data2 <- data2 %>%
  arrange(PAF_mort) %>%  # 对每个 sex 组内按 PAF_mort 排序
  mutate(site = factor(site, levels = unique(site))) # 确保排序后不改变 site 的因子级别

ggplot(data2, aes(x = PAF_mort, y = site)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, fill = "#FF99CC") +  # 设置条形颜色为绿色
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100), labels = scales::percent_format(scale = 1), expand = c(0, 0)) +  # 扩展x轴范围
  scale_y_discrete() +  # 使用 discrete scale for y axis
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.text.y = element_text(size = 9, family = "Arial", color = "black"),  # 纵坐标标签样式
    axis.line.x = element_line(color = "grey", size = 0.3),
    axis.line.y = element_line(color = "grey", size = 0.3),
    axis.ticks = element_line(color = "grey", size = 0.3),
    legend.position = "bottom",  # 将图例放置在底部
    legend.title = element_blank(),  # 删除图例标题
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_text(size = 9, family = "Arial", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
    strip.placement = "left",  # 将facet标签放置在左侧
    plot.margin = margin(l = 50),  # 调整左边距以腾出空间给文本标签
    axis.text.x = element_text(family = "Arial", color = "black")  # 设置横坐标轴字体
  ) +
  # 添加inc_attribute文本标签到x = 115的位置
  geom_text(aes(x = 115 , y = site, label = mort_95CI), 
            hjust = 0, size = 2.5, family = "Arial", color = "black") +  # 将标签放到x = 115位置
  coord_cartesian(clip = "off") +
  # 在条形图右侧添加总和标签
  geom_text(data = data2, aes(x = PAF_mort + 1, y = site, label = PAF_mort1), 
            vjust = 0.5, hjust = 0, size = 2.5, family = "Arial", color = "black") +
  # 添加"Attributable Cases"标签
  geom_text(data = data2, 
            aes(x = 150, y = 24.5, label = "Attributable Cases"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Arial", color = "black") + 
  # 添加"PAF(%)"标签
  geom_text(data = data2, 
            aes(x = 60, y = 24.5, label = "PAF(%)"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Arial", color = "black") +
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100),labels = c("0", "50", "100"), expand = c(0, 0))  # 扩展x轴范围并设置刻度标签



ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按癌种区分的PAF条图/inc-kind0-sex0.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 7.0, 
       height = 4.9, 
       units = "in", 
       dpi = 1000)



#==================================最新：按癌种区分的所有危险因素PAF条图(male+female)====================================
library(gridExtra)
library(ggplot2)
library(grid)
library(readxl)
library(dplyr)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_case_all_exposure_ICD10(all_age).xlsx")

data1 <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/所有危险因素-所有癌症PAF及case.xlsx")
data1 <- data1 %>%
  mutate(site = "All")
data1 <- data1 %>%
  mutate(ICD10 = "All")

data<-rbind(data,data1)

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 0)

data$PAF_inc<-data$PAF_inc*100
data$PAF_lower_inc<-data$PAF_lower_inc*100
data$PAF_upper_inc<-data$PAF_upper_inc*100
data$PAF_inc1 <- sprintf("%0.1f", data$PAF_inc)

data$PAF_mort<-data$PAF_mort*100
data$PAF_lower_mort<-data$PAF_lower_mort*100
data$PAF_upper_mort<-data$PAF_upper_mort*100
data$PAF_mort1 <- sprintf("%0.1f", data$PAF_mort)

data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

#inc
data1 <- data %>% filter(sex == "Male")
data1 <- data1 %>%
  arrange(PAF_inc) %>%  # 对每个 sex 组内按 PAF_inc 排序
  mutate(site = factor(site, levels = unique(site))) # 确保排序后不改变 site 的因子级别
#"Male" = "#74a9c5""#3366CC", "Female" = "#dd7389""#FF99CC" , "Total" = "#D0D4D7""#f6cd96"
plot1 <- ggplot(data1, aes(x = PAF_inc, y = site)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, fill = "#3366CC") +  # 设置条形颜色为绿色
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100), labels = scales::percent_format(scale = 1), expand = c(0, 0)) +  # 扩展x轴范围
  scale_y_discrete() +  # 使用 discrete scale for y axis
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),  # 去除x轴标签
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.text.y = element_text(size = 9, family = "Arial", color = "black"),  # 纵坐标标签样式
    axis.line.x = element_line(color = "grey", size = 0.3),
    axis.line.y = element_line(color = "grey", size = 0.3),
    axis.ticks = element_line(color = "grey", size = 0.3),
    legend.position = "bottom",  # 将图例放置在底部
    legend.title = element_blank(),  # 删除图例标题
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_text(size = 9, family = "Arial", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
    strip.placement = "left",  # 将facet标签放置在左侧
    plot.margin = margin(l = 50)  # 调整左边距以腾出空间给文本标签
  ) +
  # 添加inc_attribute文本标签到x = 115的位置
  geom_text(aes(x = 115 , y = site, label = inc_95CI), 
            hjust = 0, size = 2.5, family = "Arial", color = "black") +  # 将标签放到x = 115位置
  coord_cartesian(clip = "off") +
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = PAF_inc + 1, y = site, label = PAF_inc1), 
            vjust = 0.5, hjust = 0, size = 2.5, family = "Arial", color = "black") +
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Male", gp = gpar(fontsize = 10, fontfamily = "Arial", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -120, xmax = 0, ymin = 0, ymax = length(unique(data1$site))  # 设置文本的坐标
  )+ 
  # 只在'Male'面板上添加"Attributable Cases"标签
  geom_text(data = data1, 
            aes(x = 150, y = 20, label = "Attributable Cases"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Arial", color = "black")+ 
  # 只在'Male'面板上添加"PAF(%)"标签
  geom_text(data = data1, 
            aes(x = 60, y = 20, label = "PAF(%)"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Arial", color = "black")


data3 <- data %>% filter(sex == "Female")
data3 <- data3 %>%
  arrange(PAF_inc) %>%  # 对每个 sex 组内按 PAF_inc 排序
  mutate(site = factor(site, levels = unique(site))) # 确保排序后不改变 site 的因子级别

plot3 <- ggplot(data3, aes(x = PAF_inc, y = site)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, fill = "#FF99CC" ) +  # 设置条形颜色为灰色
  scale_x_continuous(
    limits = c(0, 150), 
    breaks = c(0, 50, 100), 
    labels = c("0", "50", "100"),  # 设置横坐标轴刻度标签
    expand = c(0, 0)
  ) +  # 扩展x轴范围并设置刻度标签
  scale_y_discrete() +  # 使用 discrete scale for y axis
  labs(x = NULL, y = NULL) +  # 删除x轴标题
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 9, family = "Arial", color = "black"),  # 横坐标标签样式
    axis.text.y = element_text(size = 9, family = "Arial", color = "black"),  # 纵坐标标签样式
    axis.line.x = element_line(color = "grey", size = 0.3),
    axis.line.y = element_line(color = "grey", size = 0.3),
    axis.ticks = element_line(color = "grey", size = 0.3),
    legend.position = "bottom",  # 将图例放置在底部
    legend.title = element_blank(),  # 删除图例标题
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_text(size = 9, family = "Arial", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
    strip.placement = "left",  # 将facet标签放置在左侧
    plot.margin = margin(l = 50)  # 调整左边距以腾出空间给文本标签
  ) +
  # 添加inc_attribute文本标签到x = 115的位置
  geom_text(aes(x = 115 , y = site, label = inc_95CI), 
            hjust = 0, size = 2.5, family = "Arial", color = "black") +  # 将标签放到x = 115位置
  coord_cartesian(clip = "off") +
  # 在条形图右侧添加总和标签
  geom_text(data = data3, aes(x = PAF_inc + 1, y = site, label = PAF_inc1), 
            vjust = 0.5, hjust = 0, size = 2.5, family = "Arial", color = "black") +
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 10, fontfamily = "Arial", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -120, xmax = 0, ymin = 0, ymax = length(unique(data3$site))  # 设置文本的坐标
  )


plot1 <- plot1 + theme(plot.margin = margin(t = 10, r = 10, b = 1, l = 50))  # Adjust for plot1
plot3 <- plot3 + theme(plot.margin = margin(t = 5, r = 10, b = 10, l = 50))  # Increase bottom margin for plot3

# 使用 grid.arrange() 将三个图排列
plot <- grid.arrange(plot1, plot3, ncol = 1, heights = c(0.9,1.1))



ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按癌种区分的PAF条图/inc-kind0-sex12.tiff", 
       plot = plot, 
       device = "tiff", 
       width = 7.1, 
       height = 5.5, 
       units = "in", 
       dpi = 1000)



#mort
data1 <- data %>% filter(sex == "Male")
data1 <- data1 %>%
  arrange(PAF_mort) %>%  # 对每个 sex 组内按 PAF_mort 排序
  mutate(site = factor(site, levels = unique(site))) # 确保排序后不改变 site 的因子级别
#"Male" = "#74a9c5""#3366CC", "Female" = "#dd7389""#FF99CC" , "Total" = "#D0D4D7""#f6cd96"
plot1 <- ggplot(data1, aes(x = PAF_mort, y = site)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, fill = "#3366CC") +  # 设置条形颜色为绿色
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100), labels = scales::percent_format(scale = 1), expand = c(0, 0)) +  # 扩展x轴范围
  scale_y_discrete() +  # 使用 discrete scale for y axis
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),  # 去除x轴标签
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.text.y = element_text(size = 9, family = "Arial", color = "black"),  # 纵坐标标签样式
    axis.line.x = element_line(color = "grey", size = 0.3),
    axis.line.y = element_line(color = "grey", size = 0.3),
    axis.ticks = element_line(color = "grey", size = 0.3),
    legend.position = "bottom",  # 将图例放置在底部
    legend.title = element_blank(),  # 删除图例标题
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_text(size = 9, family = "Arial", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
    strip.placement = "left",  # 将facet标签放置在左侧
    plot.margin = margin(l = 50)  # 调整左边距以腾出空间给文本标签
  ) +
  # 添加inc_attribute文本标签到x = 115的位置
  geom_text(aes(x = 115 , y = site, label = mort_95CI), 
            hjust = 0, size = 2.5, family = "Arial", color = "black") +  # 将标签放到x = 115位置
  coord_cartesian(clip = "off") +
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = PAF_mort + 1, y = site, label = PAF_mort1), 
            vjust = 0.5, hjust = 0, size = 2.5, family = "Arial", color = "black") +
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Male", gp = gpar(fontsize = 10, fontfamily = "Arial", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -120, xmax = 0, ymin = 0, ymax = length(unique(data1$site))  # 设置文本的坐标
  )+ 
  # 只在'Male'面板上添加"Attributable Cases"标签
  geom_text(data = data1, 
            aes(x = 150, y = 20, label = "Attributable Cases"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Arial", color = "black")+ 
  # 只在'Male'面板上添加"PAF(%)"标签
  geom_text(data = data1, 
            aes(x = 60, y = 20, label = "PAF(%)"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Arial", color = "black")


data3 <- data %>% filter(sex == "Female")
data3 <- data3 %>%
  arrange(PAF_mort) %>%  # 对每个 sex 组内按 PAF_mort 排序
  mutate(site = factor(site, levels = unique(site))) # 确保排序后不改变 site 的因子级别

plot3 <- ggplot(data3, aes(x = PAF_mort, y = site)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, fill = "#FF99CC" ) +  # 设置条形颜色为灰色
  scale_x_continuous(
    limits = c(0, 150), 
    breaks = c(0, 50, 100), 
    labels = c("0", "50", "100"),  # 设置横坐标轴刻度标签
    expand = c(0, 0)
  ) +  # 扩展x轴范围并设置刻度标签
  scale_y_discrete() +  # 使用 discrete scale for y axis
  labs(x = NULL, y = NULL) +  # 删除x轴标题
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 9, family = "Arial", color = "black"),  # 横坐标标签样式
    axis.text.y = element_text(size = 9, family = "Arial", color = "black"),  # 纵坐标标签样式
    axis.line.x = element_line(color = "grey", size = 0.3),
    axis.line.y = element_line(color = "grey", size = 0.3),
    axis.ticks = element_line(color = "grey", size = 0.3),
    legend.position = "bottom",  # 将图例放置在底部
    legend.title = element_blank(),  # 删除图例标题
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_text(size = 9, family = "Arial", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
    strip.placement = "left",  # 将facet标签放置在左侧
    plot.margin = margin(l = 50)  # 调整左边距以腾出空间给文本标签
  ) +
  # 添加inc_attribute文本标签到x = 115的位置
  geom_text(aes(x = 115 , y = site, label = mort_95CI), 
            hjust = 0, size = 2.5, family = "Arial", color = "black") +  # 将标签放到x = 115位置
  coord_cartesian(clip = "off") +
  # 在条形图右侧添加总和标签
  geom_text(data = data3, aes(x = PAF_mort + 1, y = site, label = PAF_mort1), 
            vjust = 0.5, hjust = 0, size = 2.5, family = "Arial", color = "black") +
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 10, fontfamily = "Arial", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -120, xmax = 0, ymin = 0, ymax = length(unique(data3$site))  # 设置文本的坐标
  )


plot1 <- plot1 + theme(plot.margin = margin(t = 10, r = 10, b = 1, l = 50))  # Adjust for plot1
plot3 <- plot3 + theme(plot.margin = margin(t = 5, r = 10, b = 10, l = 50))  # Increase bottom margin for plot3

# 使用 grid.arrange() 将三个图排列
plot <- grid.arrange(plot1, plot3, ncol = 1, heights = c(0.9,1.1))



ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按癌种区分的PAF条图/mort-kind0-sex12.tiff", 
       plot = plot, 
       device = "tiff", 
       width = 7.1, 
       height = 5.5, 
       units = "in", 
       dpi = 1000)





#============================按exposure区分的各癌种case条图inc(仅有男女合计）（更新）==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)
library(grid)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")

data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 0)

# 定义癌种的顺序
site_order <- c(
  "Traches,bronchus and lung","Nasopharynx","Oral cavity, pharynx","Larynx",
  "Liver","Colorectum","Stomach","Esophagus","Pancreas","Gallbladder",
  "Cervix uteri","Breast","Corpus uteri","Prostate","Ovary","Vulva, Vagina", "Penis",
  "Thyroid","Bladder","Kidney","Myeloid leukaemia","Hodgkin disease")

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
data <- data %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序

custom_colors <- c("#567CD0","#82A3E0","#A9C2EA" ,"#D6E8FA", 
                   "#F2A7C0", "#F5B6CC", "#F8C5D8", "#FBD4E4", "#FEE3ED", "#FFF2F7", 
                   "#E47A6D", "#E68E80", "#E8A292", "#EAB6A4", "#EDCAB6", "#EFD0A8", "#F1D6B9",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8")


data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

data <- data %>%
  mutate(sex = as.character(sex))

data <- data %>%
  group_by(exposure, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(inc_attribute))

data1 <- data %>% filter(sex == "Total")

exposure_order1 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Low vegetable intake","Helicobacter pylori","Second-hand smoking",
  "Red meat","Excess bodyweight","Human papillomavirus","Epstein-Barr virus",
  "Hepatitis B virus", "Low fruit intake","Diabetes", "Alcohol","Smoking")
data1 <- data1 %>% mutate(exposure = factor(exposure, levels = exposure_order1))  # 按照提供的顺序调整 site 的因子顺序

ggplot(data1, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 26000), breaks = c(0,5000, 10000, 15000, 20000, 25000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = NULL, y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Arial", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_text(size = 9, family = "Arial", face = "bold", color = "black"),  # 设置图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Arial", color = "black"),  # 图例文本样式
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示



ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按exposure区分的PAF条图/inc-kind0-sex0.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 12.7, 
       height = 4.2, 
       units = "in", 
       dpi = 1000)

#============================按exposure区分的各癌种case条图mort(仅有男女合计）（更新）==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)
library(grid)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)

data <- data %>% filter(kind == 0)

# 定义癌种的顺序
exposure_order1 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Low vegetable intake","Helicobacter pylori","Second-hand smoking",
  "Human papillomavirus","Excess bodyweight","Red meat","Epstein-Barr virus",
  "Diabetes", "Low fruit intake", "Alcohol","Hepatitis B virus","Smoking")

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
data <- data %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序

custom_colors <- c("#567CD0","#82A3E0","#A9C2EA" ,"#D6E8FA", 
                   "#F2A7C0", "#F5B6CC", "#F8C5D8", "#FBD4E4", "#FEE3ED", "#FFF2F7", 
                   "#E47A6D", "#E68E80", "#E8A292", "#EAB6A4", "#EDCAB6", "#EFD0A8", "#F1D6B9",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8")


data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

data <- data %>%
  mutate(sex = as.character(sex))

data <- data %>%
  group_by(exposure, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(mort_attribute))

data1 <- data %>% filter(sex == "Total")

exposure_order1 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Low vegetable intake","Helicobacter pylori","Second-hand smoking",
  "Human papillomavirus","Excess bodyweight","Red meat","Epstein-Barr virus",
  "Diabetes", "Low fruit intake", "Alcohol","Hepatitis B virus","Smoking")
data1 <- data1 %>% mutate(exposure = factor(exposure, levels = exposure_order1))  # 按照提供的顺序调整 site 的因子顺序

ggplot(data1, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 20000), breaks = c(0,3000, 6000, 9000, 12000, 15000,18000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = NULL, y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Arial", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Arial", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示


ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按exposure区分的PAF条图/mort-kind0-sex0.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 12.7, 
       height = 4.2, 
       units = "in", 
       dpi = 1000)


#============================按exposure区分的各癌种case条图inc(male+female)（更新）==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)
library(grid)
library(gridExtra)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 0)

# 定义癌种的顺序
site_order <- c("Traches,bronchus and lung","Nasopharynx","Oral cavity, pharynx","Larynx",
                "Liver","Colorectum","Stomach","Esophagus","Pancreas","Gallbladder",
                "Cervix uteri","Breast","Corpus uteri","Prostate","Ovary","Vulva, Vagina", "Penis",
                "Thyroid","Bladder","Kidney","Myeloid leukaemia","Hodgkin disease")

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
data <- data %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序


custom_colors <- c("#567CD0","#82A3E0","#A9C2EA" ,"#D6E8FA", 
                   "#F2A7C0", "#F5B6CC", "#F8C5D8", "#FBD4E4", "#FEE3ED", "#FFF2F7", 
                   "#E47A6D", "#E68E80", "#E8A292", "#EAB6A4", "#EDCAB6", "#EFD0A8", "#F1D6B9",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8")


data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

data <- data %>%
  mutate(sex = as.character(sex))

data <- data %>%
  group_by(exposure, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(inc_attribute))  

plot0 <- ggplot(data, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 26000), breaks = c(0,5000, 10000, 15000, 20000, 25000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Arial", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Arial", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_text(size = 9, family = "Arial", face = "bold", color = "black"),  # 设置图例标题
        legend.text = element_text(size = 8, family = "Arial", color = "black"),  # 设置图例文本样式
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Arial", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -12000, xmax = -10000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示


data1 <- data %>% filter(sex == "Male")
exposure_order1 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Low vegetable intake","Second-hand smoking",
  "Human papillomavirus","Helicobacter pylori","Excess bodyweight","Red meat","Epstein-Barr virus",
  "Diabetes", "Low fruit intake","Hepatitis B virus", "Alcohol","Smoking"
)
data1 <- data1 %>%
  mutate(exposure = factor(exposure, levels = exposure_order1))  # 按照提供的顺序调整 site 的因子顺序

plot1 <- ggplot(data1, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 26000), breaks = c(0,5000, 10000, 15000, 20000, 25000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = NULL, y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 去除x轴标签
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.text.y = element_text(size = 9, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Arial", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Male", gp = gpar(fontsize = 11, fontfamily = "Arial", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -6000, xmax = -2800, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示


data3 <- data %>% filter(sex == "Female")
exposure_order3 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Low vegetable intake","Physical inactivity","Helicobacter pylori",
  "Smoking","Red meat","Second-hand smoking","Hepatitis B virus","Alcohol","Epstein-Barr virus",
  "Excess bodyweight","Low fruit intake","Diabetes","Human papillomavirus"
)
data3 <- data3 %>%
  mutate(exposure = factor(exposure, levels = exposure_order3))  # 按照提供的顺序调整 site 的因子顺序

plot3 <- ggplot(data3, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 26000), breaks = c(0,5000, 10000, 15000, 20000, 25000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Arial", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Arial", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.text = element_blank(),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Arial", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -6000, xmax = -2800, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data3, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

plot1 <- plot1 + theme(plot.margin = margin(t = 10, r = 100, b = 1, l = 50))  # Adjust for plot1
plot3 <- plot3 + theme(plot.margin = margin(t = 5, r = 100, b = 1, l = 50))  # Increase bottom margin for plot3

# 使用 grid.arrange() 将三个图排列
plot <- grid.arrange(plot1, plot3, ncol = 1, heights = c(1,1))

library(cowplot)

# 提取plot2的图例
legend <- get_legend(plot0)

# 将图例添加到右侧
final_plot <- plot_grid(plot, legend, ncol = 2, rel_widths = c(0.9, 0.1))

# 显示最终的组合图
final_plot


ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按exposure区分的PAF条图/inc_sex12_kind0.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 15,
       height = 6, 
       units = "in", 
       dpi = 1000)



#============================按exposure区分的各癌种case条图mort(male+female)（更新）==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)
library(grid)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)

data <- data %>% filter(kind == 0)

# 定义癌种的顺序
site_order <- c("Traches,bronchus and lung","Nasopharynx","Oral cavity, pharynx","Larynx",
                "Liver","Colorectum","Stomach","Esophagus","Pancreas","Gallbladder",
                "Cervix uteri","Breast","Corpus uteri","Prostate","Ovary","Vulva, Vagina", "Penis",
                "Thyroid","Bladder","Kidney","Myeloid leukaemia","Hodgkin disease")

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
data <- data %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序


custom_colors <- c("#567CD0","#82A3E0","#A9C2EA" ,"#D6E8FA", 
                   "#F2A7C0", "#F5B6CC", "#F8C5D8", "#FBD4E4", "#FEE3ED", "#FFF2F7", 
                   "#E47A6D", "#E68E80", "#E8A292", "#EAB6A4", "#EDCAB6", "#EFD0A8", "#F1D6B9",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8")


data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

data <- data %>%
  mutate(sex = as.character(sex))

data <- data %>%
  group_by(exposure, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(mort_attribute))  

plot0 <- ggplot(data, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 26000), breaks = c(0,5000, 10000, 15000, 20000, 25000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Arial", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Arial", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_text(size = 9, family = "Arial", face = "bold", color = "black"),  # 设置图例标题
        legend.text = element_text(size = 8, family = "Arial", color = "black"),  # 设置图例文本样式
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Arial", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -12000, xmax = -10000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示


data1 <- data %>% filter(sex == "Male")
exposure_order1 <- c(
  "Physical inactivity","Clonorchis sinensis","Hepatitis C virus","Low vegetable intake",
  "Human papillomavirus","Helicobacter pylori","Second-hand smoking","Excess bodyweight","Red meat","Epstein-Barr virus",
  "Diabetes", "Low fruit intake","Hepatitis B virus", "Alcohol","Smoking"
)

data1 <- data1 %>%
  mutate(exposure = factor(exposure, levels = exposure_order1))  # 按照提供的顺序调整 site 的因子顺序

plot1 <- ggplot(data1, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 26000), breaks = c(0,5000, 10000, 15000, 20000, 25000),expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = NULL, y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 去除x轴标签
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.text.y = element_text(size = 9, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Arial", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Male", gp = gpar(fontsize = 11, fontfamily = "Arial", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -6000, xmax = -2800, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示


data3 <- data %>% filter(sex == "Female")
exposure_order3 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Low vegetable intake","Physical inactivity","Helicobacter pylori",
  "Smoking","Red meat","Second-hand smoking","Alcohol","Excess bodyweight","Epstein-Barr virus",
  "Hepatitis B virus","Diabetes","Human papillomavirus","Low fruit intake"
)
data3 <- data3 %>%
  mutate(exposure = factor(exposure, levels = exposure_order3))  # 按照提供的顺序调整 site 的因子顺序

plot3 <- ggplot(data3, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 26000), breaks = c(0,5000, 10000, 15000, 20000, 25000),expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Arial", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Arial", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.text = element_blank(),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Arial", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -6000, xmax = -2800, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data3, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

plot1 <- plot1 + theme(plot.margin = margin(t = 10, r = 100, b = 1, l = 50))  # Adjust for plot1
plot3 <- plot3 + theme(plot.margin = margin(t = 5, r = 100, b = 1, l = 50))  # mortrease bottom margin for plot3

# 使用 grid.arrange() 将三个图排列
plot <- grid.arrange(plot1, plot3, ncol = 1, heights = c(1,1))

library(cowplot)

# 提取plot2的图例
legend <- get_legend(plot0)

# 将图例添加到右侧
final_plot <- plot_grid(plot, legend, ncol = 2, rel_widths = c(0.9, 0.1))

# 显示最终的组合图
final_plot


ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按exposure区分的PAF条图/mort_sex12_kind0.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 15, 
       height = 6, 
       units = "in", 
       dpi = 1000)



#===================最新-按factor区分的PAF条图(total)==============================================
library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)
library(scales)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/某一factor所有癌症PAF及case.xlsx")
data1 <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/所有危险因素-所有癌症PAF及case.xlsx")
data1 <- data1 %>%
  mutate(factor = "All")
data<-rbind(data,data1)

data <- data %>% filter(kind == 0)

data$PAF_inc<-data$PAF_inc*100
data$PAF_mort<-data$PAF_mort*100

# 定义暴露因素的顺序
factor_order <- c("Metabolic factor","Dietary factor","Infectious agent","Behaviour factor","All")

# 确保每个癌种都包含所有暴露因素，并按暴露因素的顺序排列
data <- data %>%
  mutate(factor = factor(factor, levels = factor_order)
  )  # 确保暴露因子按顺序排列

data$sex <- factor(data$sex, levels = c(2, 1, 0), labels = c("Female", "Male","Total" ))

#inc
ggplot(data, aes(y = factor, x = PAF_inc, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +  # 绘制条形图
  scale_x_continuous(name = "PAF(%)", 
                     limits = c(0, 60),
                     breaks = c(0, 20, 40, 60), 
                     labels = scales::comma_format(),  # 去掉百分号，改为普通数字格式
                     expand = c(0, 0)) +  # 横坐标范围
  scale_fill_manual(values = c("Male" = "#3366CC", "Female" = "#FF99CC" , "Total" = "#D0D4D7")) +
  labs(x = "PAF(%)", y = NULL) + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 10, family = "Arial", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 10, family = "Arial", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.position = "right",  # 图例放置在右侧
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 9, family = "Arial", color = "black"),  # 图例文本样式
        panel.grid = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside") +  # 设置面板标签在纵坐标轴外侧
  geom_text(aes(x = PAF_inc + 1, y = factor, label = sprintf("%.1f", PAF_inc)),  # 保留一位小数
            position = position_dodge(width = 0.8),  # 保证标签的位置与条形图对齐
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black") +
  # 在x = 100的位置添加inc_95CI文本标签
  geom_text(aes(x = 105, y = factor, label = inc_95CI),  # 设置x = 100，确保文本标签在横坐标轴外侧
            position = position_dodge(width = 0.8),  # 确保文本与条形图对齐
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black") 

ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按factor区分的PAF条图/inc-Both_kind.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 6, 
       height = 3, 
       units = "in", 
       dpi = 1000)

#mort
ggplot(data, aes(y = factor, x = PAF_mort, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +  # 绘制条形图
  scale_x_continuous(name = "PAF(%)", 
                     limits = c(0, 60),
                     breaks = c(0, 20, 40, 60), 
                     labels = scales::comma_format(),  # 去掉百分号，改为普通数字格式
                     expand = c(0, 0)) +  # 横坐标范围
  scale_fill_manual(values = c("Male" = "#3366CC", "Female" = "#FF99CC" , "Total" = "#D0D4D7")) +
  labs(x = "PAF(%)", y = NULL) + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 10, family = "Arial", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 10, family = "Arial", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.position = "right",  # 图例放置在右侧
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 9, family = "Arial", color = "black"),  # 图例文本样式
        panel.grid = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside") +  # 设置面板标签在纵坐标轴外侧
  geom_text(aes(x = PAF_mort + 1, y = factor, label = sprintf("%.1f", PAF_mort)),  # 保留一位小数
            position = position_dodge(width = 0.8),  # 保证标签的位置与条形图对齐
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black") +
  # 在x = 100的位置添加mort_95CI文本标签
  geom_text(aes(x = 105, y = factor, label = mort_95CI),  # 设置x = 100，确保文本标签在横坐标轴外侧
            position = position_dodge(width = 0.8),  # 确保文本与条形图对齐
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black") 

ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按factor区分的PAF条图/mort-Both_kind.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 6, 
       height = 3, 
       units = "in", 
       dpi = 1000)

#===================最新-按factor区分的PAF条图(合并Urban、Rural)==============================================
library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)
library(scales)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/某一factor所有癌症PAF及case.xlsx")
data1 <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/所有危险因素-所有癌症PAF及case.xlsx")
data1 <- data1 %>%
  mutate(factor = "All")
data<-rbind(data,data1)

data <- data %>% filter(kind != 0)

data$PAF_inc<-data$PAF_inc*100
data$PAF_mort<-data$PAF_mort*100

# 定义暴露因素的顺序
factor_order <- c("All","Infectious agent", "Metabolic factor","Dietary factor","Behaviour factor")

# 确保每个癌种都包含所有暴露因素，并按暴露因素的顺序排列
data <- data %>%
  mutate(factor = factor(factor, levels = factor_order)
  )  # 确保暴露因子按顺序排列

data$sex <- factor(data$sex, levels = c(2, 1, 0), labels = c("Female", "Male","Total" ))

#"Male" = "#7fb2d5", "Female" = "#f47f72" , "Total" = "#f6cd96""#fbd178""#abd9d0""#ABADAF""#f6cd96"

panel_labels <- c("1" = "Urban", "2" = "Rural")

#inc
ggplot(data, aes(y = factor, x = PAF_inc, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +  # 绘制条形图
  scale_x_continuous(name = "PAF(%)", 
                     limits = c(0, 60),
                     breaks = c(0, 20, 40, 60), 
                     labels = scales::comma_format(),  # 去掉百分号，改为普通数字格式
                     expand = c(0, 0)) +  # 横坐标范围
  scale_fill_manual(values = c("Male" = "#3366CC", "Female" = "#FF99CC" , "Total" = "#D0D4D7")) +
  labs(x = "PAF(%)", y = NULL) + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 10, family = "Arial", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 10, family = "Arial", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.position = "right",  # 图例放置在右侧
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 9, family = "Arial", color = "black"),  # 图例文本样式
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside") +  # 设置面板标签在纵坐标轴外侧
  geom_text(aes(x = PAF_inc + 1, y = factor, label = sprintf("%.1f", PAF_inc)),  # 保留一位小数
            position = position_dodge(width = 0.8),  # 保证标签的位置与条形图对齐
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black") +
  # 在x = 100的位置添加inc_95CI文本标签
  geom_text(aes(x = 105, y = factor, label = inc_95CI),  # 设置x = 100，确保文本标签在横坐标轴外侧
            position = position_dodge(width = 0.8),  # 确保文本与条形图对齐
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black") +
  facet_wrap(~ kind, ncol = 1, scales = "free_y", 
             labeller = as_labeller(panel_labels), 
             strip.position = "left")  # 将面板标签放置在左侧

ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按factor区分的PAF条图/inc-kind12.tiff",  
       plot = last_plot(), 
       device = "tiff", 
       width = 8, 
       height = 5, 
       units = "in", 
       dpi = 1000)


#mort
ggplot(data, aes(y = factor, x = PAF_mort, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +  # 绘制条形图
  scale_x_continuous(name = "PAF(%)", 
                     limits = c(0, 60),
                     breaks = c(0, 20, 40, 60), 
                     labels = scales::comma_format(),  # 去掉百分号，改为普通数字格式
                     expand = c(0, 0)) +  # 横坐标范围
  scale_fill_manual(values = c("Male" = "#3366CC", "Female" = "#FF99CC" , "Total" = "#D0D4D7")) +
  labs(x = "PAF(%)", y = NULL) + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, family = "Arial", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 10, family = "Arial", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 10, family = "Arial", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.position = "right",  # 图例放置在右侧
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 9, family = "Arial", color = "black"),  # 图例文本样式
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Arial", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside") +  # 设置面板标签在纵坐标轴外侧
  geom_text(aes(x = PAF_mort + 1, y = factor, label = sprintf("%.1f", PAF_mort)),  # 保留一位小数
            position = position_dodge(width = 0.8),  # 保证标签的位置与条形图对齐
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black") +
  # 在x = 100的位置添加inc_95CI文本标签
  geom_text(aes(x = 105, y = factor, label = mort_95CI),  # 设置x = 100，确保文本标签在横坐标轴外侧
            position = position_dodge(width = 0.8),  # 确保文本与条形图对齐
            vjust = 0.5, hjust = 0, size = 3, family = "Arial", color = "black") +
  facet_wrap(~ kind, ncol = 1, scales = "free_y", 
             labeller = as_labeller(panel_labels), 
             strip.position = "left")  # 将面板标签放置在左侧

ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按factor区分的PAF条图/mort-kind12.tiff",  
       plot = last_plot(), 
       device = "tiff", 
       width = 8, 
       height = 5, 
       units = "in", 
       dpi = 1000)


#==========================纵坐标site+exposure，横坐标All\sex\kind\age的PAF热图==============================================

install.packages("ggtext")
install.packages("stringr")

library(ggplot2)
library(dplyr)
library(readxl)  
library(ggtext) 
library(stringr)
library(scales)


data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_age_kind012_sex012.xlsx")
data1<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")
data2<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_all_exposure_ICD10.xlsx")
data3<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_all_exposure_ICD10(all_age).xlsx")


data1 <- data1 %>%
  mutate(age = "All")
data<-rbind(data,data1)
data<-data[,-c(9:14,21:24)]

data3 <- data3 %>%
  mutate(age = "All")
data2<-rbind(data2,data3)
data2<-data2[,-c(14,15)]
data2 <- data2 %>%
  mutate(exposure = "All")

data<-rbind(data,data2)

data$PAF_inc <- as.numeric(data$PAF_inc)
data$PAF_lower_inc <- as.numeric(data$PAF_lower_inc)
data$PAF_upper_inc <- as.numeric(data$PAF_upper_inc)

data$inc_attribute<-data$inc_all*data$PAF_inc
data$inc_attribute_lower<-data$inc_all*data$PAF_lower_inc
data$inc_attribute_upper<-data$inc_all*data$PAF_upper_inc

data$PAF_mort <- as.numeric(data$PAF_mort)
data$PAF_lower_mort <- as.numeric(data$PAF_lower_mort)
data$PAF_upper_mort <- as.numeric(data$PAF_upper_mort)

data$mort_attribute<-data$mort_all*data$PAF_mort
data$mort_attribute_lower<-data$mort_all*data$PAF_lower_mort
data$mort_attribute_upper<-data$mort_all*data$PAF_upper_mort



#data <- data %>%
# mutate(site = str_replace(site, "Gallbladder etc.", "Gallbladder"))

#data$age <- ifelse(data$age == "all", "All", data$age)

data1 <- data %>% filter(kind == 0 & sex ==1 & age == "All")
data2 <- data %>% filter(kind == 0 & sex ==2 & age == "All")
data5 <- data %>% filter(kind == 0 & sex ==0)
data5 <- data5 %>%
  filter(age == "20-29" | age == "30-39")
data5<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+sex+kind+ICD10+site,data5,sum)
data6 <- data %>% filter(kind == 0 & sex ==0)
data6 <- data6 %>%
  filter(age == "60-69" | age == "70-")
data6<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+sex+kind+ICD10+site,data6,sum)
data7<-data %>% filter(kind == 0 & sex ==0 & age == "All")
data8 <- data %>% filter(kind == 0 & sex ==0)
data8 <- data8 %>%
  filter(age =="40-49" | age == "50-59")
data8<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+sex+kind+ICD10+site,data8,sum)


data1<-data1[,-c(2:4,9:14)]
data2<-data2[,-c(2:4,9:14)]
data3<-data3[,-c(2:4,9:14)]
data4<-data4[,-c(2:4,9:14)]
data7<-data7[,-c(2:4,9:14)]
data5<-data5[,-c(2:3)]
data6<-data6[,-c(2:3)]
data8<-data8[,-c(2:3)]
data1 <- data1 %>%
  mutate(factor = "Male")
data2 <- data2 %>%
  mutate(factor = "Female")
data5 <- data5 %>%
  mutate(factor = "Age=20-39")
data6 <- data6 %>%
  mutate(factor = "Age≥60")
data7 <- data7 %>%
  mutate(factor = "All")
data8 <- data8 %>%
  mutate(factor = "Age=40-59")

data<-rbind(data1,data2,data5,data6,data7,data8)

data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all


#填充age20-39的Penis
new_row <- data.frame(
  exposure = "Human papillomavirus",
  site = "Penis",
  sex = 0,
  kind = 0,
  ICD10 = "C60",
  PAF_inc = 0.442,
  PAF_mort = 0.442,
  factor = "Age=20-39",
  stringsAsFactors = FALSE  # 避免将字符向量转换为因子
)
# 添加其他列，这些列的值为 NA
other_columns <- setdiff(names(data), names(new_row))
new_row[other_columns] <- NA

# 将新行添加到原始数据框中
data <- bind_rows(data, new_row)

new_row <- data.frame(
  exposure = "All",
  site = "Penis",
  sex = 0,
  kind = 0,
  ICD10 = "C60",
  PAF_inc = 0.442,
  PAF_mort = 0.442,
  factor = "Age=20-39",
  stringsAsFactors = FALSE  # 避免将字符向量转换为因子
)
# 添加其他列，这些列的值为 NA
other_columns <- setdiff(names(data), names(new_row))
new_row[other_columns] <- NA

# 将新行添加到原始数据框中
data <- bind_rows(data, new_row)


#
data$PAF_inc<-data$PAF_inc*100
data$PAF_mort<-data$PAF_mort*100

data <- data %>%
  mutate(
    site = factor(site, levels = c(
      "Penis",  "Prostate", "Vulva, Vagina", "Breast", "Ovary", "Corpus uteri", "Cervix uteri",
      "Non-hodgkin lymphoma", "Hodgkin disease", "Myeloid leukaemia", "Thyroid", "Bladder", "Kidney",
      "Gallbladder", "Liver", "Pancreas", "Colorectum", "Stomach", "Esophagus", "Traches,bronchus and lung",
      "Larynx", "Nasopharynx", "Oral cavity, pharynx"
    )),  # 设置 site 排序
    
    exposure = factor(exposure, levels = c(
      "Clonorchis sinensis","Hepatitis C virus", "Hepatitis B virus", "Helicobacter pylori", "Human papillomavirus", 
      "Epstein-Barr virus", "Diabetes", "Excess bodyweight", "Red meat", "Low fruit intake", 
      "Low vegetable intake", "Physical inactivity", "Alcohol", "Second-hand smoking", "Smoking","All"
    ))  # 设置 exposure 排序
  )

na_rows <- data[is.na(data$PAF_inc), ]
#data <- na.omit(data)
data$PAF_inc <- ifelse(is.na(data$PAF_inc), 0, data$PAF_inc)
data$PAF_mort <- ifelse(is.na(data$PAF_mort), 0, data$PAF_mort)
colors <- c("#3B5DA0","#3F62A4","#4970AF","#5783B9","#6697C3","#A4C3DD","#CADEEC","#EBF3F8","#FFFFFF","#FDF7DA","#FBE99D","#FAD789","#F8BB6A","#E7995C","#E18D57","#C14336","#AB1F2A")
colors <- c("#567BA9", "white", "#B74A4A")  # 从蓝色到白色再到红色

# PAF_inc
# 转换为均匀分布的值，并将 0 转换为 NA
data <- data %>%
  mutate(PAF_uniform = (rank(PAF_inc) - 1) / (n() - 1),
         PAF_uniform = ifelse(PAF_inc == 0, NA, PAF_uniform))

ggplot(data, aes(x = factor, y = interaction(exposure, site, drop = TRUE), fill = ifelse(is.na(PAF_uniform), NA, PAF_uniform))) +
  geom_tile(color = "white", width = 1, height = 1) +  # 绘制热图格子，格子之间有白色边框
  geom_text(aes(label = ifelse(!is.na(PAF_inc) & PAF_inc != 0, sprintf("%.1f", PAF_inc), "")), family = "Times New Roman", color = "black", size = 2.5) +  # 在每个格子中加入PAF值
  scale_fill_gradientn(
    colors = colors,  # 使用定义的颜色数组
    values = scales::rescale(seq(0, 1, length.out = 100)),  # 使用100个色阶进行过渡
    na.value = "grey",  # 当值为 NA 时填充灰色
    guide = guide_colorbar(
      barwidth = 8,  # 图例色条的宽度
      barheight = 0.4,  # 图例色条的高度
      ticks = FALSE,  # 不显示刻度线
      label.position = "bottom",  # 标签放在色条下方
      labels = percent_format(scale = 1)(rescale(range(data$PAF_inc, na.rm = TRUE))),  # 使用 PAF_inc 的实际值范围
      direction = "horizontal"  # 设置图例为水平排列
    )
  ) +
  labs(x = NULL, y = NULL, fill = "PAF") +  # 添加图例标题
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 9, family = "Arial", color = "black"),
    axis.line.x = element_line(color = "grey", linewidth = 0.3),
    axis.line.y = element_line(color = "grey", linewidth = 0.3),
    axis.title.x = element_text(size = 9, family = "Arial", face = "bold", color = "black"),
    axis.ticks = element_line(color = "grey", linewidth = 0.3),
    legend.title = element_text(size = 7, family = "Arial", face = "bold"),  # 设置图例标题样式
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.background = element_rect(color = NA),
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_blank(),  # 隐藏图例文本标签
    panel.grid = element_blank(),  # 去除背景网格
    strip.text = element_text(size = 5),  # 调整分面标题的字体大小
    legend.position = "bottom",  # 将图例位置设置为底部
    panel.background = element_rect(fill = "grey90"),  # 设置坐标轴内背景为灰色
    plot.background = element_rect(fill = "white", color = NA),  # 去除整个图形的边框
    plot.margin = margin(0, 0, 0, 0)  # 移除图形的外边距
  ) +
  scale_x_discrete(limits = c("All", "Male", "Female", "Age=20-39", "Age=40-59", "Age≥60"), expand = c(0, 0)) +
  scale_y_discrete(
    labels = function(x) {
      sapply(x, function(label) {
        parts <- unlist(strsplit(label, "\\."))
        if (!is.na(parts[1]) && parts[1] == "All") {
          return(paste0("<b>", parts[2], "</b>"))  # 将 site 部分加粗
        } else {
          return(parts[1])
        }
      })
    }
  ) +
  theme(
    axis.text.y = element_markdown()  # 使用 ggtext 来支持 HTML 标签
  )

# PAF_mort
data <- data %>%
  mutate(PAF_uniform = (rank(PAF_mort) - 1) / (n() - 1),
         PAF_uniform = ifelse(PAF_mort == 0, NA, PAF_uniform))

ggplot(data, aes(x = factor, y = interaction(exposure, site, drop = TRUE), fill = ifelse(is.na(PAF_uniform), NA, PAF_uniform))) +
  geom_tile(color = "white", width = 1, height = 1) +  # 绘制热图格子，格子之间有白色边框
  geom_text(aes(label = ifelse(!is.na(PAF_mort) & PAF_inc != 0, sprintf("%.1f", PAF_mort), "")), family = "Times New Roman", color = "black", size = 2.5) +  # 在每个格子中加入PAF值
  scale_fill_gradientn(
    colors = colors,  # 使用定义的颜色数组
    values = rescale(seq(0, 1, length.out = length(colors))),  # 均匀分布的值
    na.value = "grey",  # 当值为 NA 时填充灰色
    guide = guide_colorbar(
      barwidth = 8,  # 图例色条的宽度
      barheight = 0.4,  # 图例色条的高度
      ticks = FALSE,  # 不显示刻度线
      label.position = "bottom",  # 标签放在色条下方
      labels = percent_format(scale = 1)(rescale(range(data$PAF_mort, na.rm = TRUE))),  # 使用 PAF_inc 的实际值范围
      direction = "horizontal"  # 设置图例为水平排列
    )
  ) +
  labs(x = NULL, y = NULL, fill = "PAF") +  # 添加图例标题
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 9, family = "Arial", color = "black"),
    axis.line.x = element_line(color = "grey", linewidth = 0.3),
    axis.line.y = element_line(color = "grey", linewidth = 0.3),
    axis.title.x = element_text(size = 9, family = "Arial", face = "bold", color = "black"),
    axis.ticks = element_line(color = "grey", linewidth = 0.3),
    legend.title = element_text(size = 7, family = "Arial", face = "bold"),  # 设置图例标题样式
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.background = element_rect(color = NA),
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_blank(),  # 隐藏图例文本标签
    panel.grid = element_blank(),  # 去除背景网格
    strip.text = element_text(size = 5),  # 调整分面标题的字体大小
    legend.position = "bottom",  # 将图例位置设置为底部
    panel.background = element_rect(fill = "grey90"),  # 设置坐标轴内背景为灰色
    plot.background = element_rect(fill = "white", color = NA),  # 去除整个图形的边框
    plot.margin = margin(0, 0, 0, 0)  # 移除图形的外边距
  ) +
  scale_x_discrete(limits = c("All", "Male", "Female", "Age=20-39", "Age=40-59", "Age≥60"), expand = c(0, 0)) +
  scale_y_discrete(
    labels = function(x) {
      sapply(x, function(label) {
        parts <- unlist(strsplit(label, "\\."))
        if (!is.na(parts[1]) && parts[1] == "All") {
          return(paste0("<b>", parts[2], "</b>"))  # 将 site 部分加粗
        } else {
          return(parts[1])
        }
      })
    }
  ) +
  theme(
    axis.text.y = element_markdown()  # 使用 ggtext 来支持 HTML 标签
  )

ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/热图/mort.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 3.5, 
       height = 14, 
       units = "in", 
       dpi = 1200)





data1 <- read_excel("吉林3776+江苏4062+天津3837+北京1844-无重复ID.xlsx")
data2 <- read_excel("【】‘/微量元素/检测结果/汇总/吉林2892+江苏3370+天津3198-硒.xlsx")

data1 <- data1[-c(1),]
data2 <- data2[-c(1),]

names(data1)[names(data1) == "样品编号"] <- "id"
names(data2)[names(data2) == "样品编号"] <- "id"


data1_filtered <- data1 %>% 
  filter(!(id %in% data2$id))

data1_filtered <- data1_filtered[,-c(1)]


write.xlsx(data1_filtered, "~/Desktop/吉林884+江苏692+天津639+北京1844-合计4059.xlsx")


#
data <- read_excel("吉林3752+江苏3968+天津3711+北京1076-无重复ID.xlsx")

data <- data[-c(1),]
data <- data[,-c(1,2)]

names(data)[names(data) == "exposure"] <- "risk_factor"



data[, 2:26] <- lapply(data[, 2:26], function(x) as.numeric(as.character(x)))

calculate_percentiles <- function(column) {
  p2_5 <- quantile(column, 0.025, na.rm = TRUE)  
  p25 <- quantile(column, 0.25, na.rm = TRUE) 
  median_value <- median(column, na.rm = TRUE)
  p75 <- quantile(column, 0.75, na.rm = TRUE) 
  p97_5 <- quantile(column, 0.975, na.rm = TRUE)
  return(c(p2_5, p25, median_value, p75, p97_5))
}

calculate_stats <- function(x) {
  x <- as.numeric(x)
  n <- sum(!is.na(x))                 # 有效样本数
  mean_x <- mean(x, na.rm = TRUE)     # 平均值
  sd_x <- sd(x, na.rm = TRUE)         # 标准差
  se <- sd_x / sqrt(n)                # 标准误差
  ci_lower <- mean_x - 1.96 * se      # 95%CI下限
  ci_upper <- mean_x + 1.96 * se      # 95%CI上限
  return(c(Mean = mean_x, SD = sd_x, SE = se, Lower95CI = ci_lower, Upper95CI = ci_upper))
}


selected_columns <- data[, c(2:26)]
results <- apply(selected_columns, 2, calculate_stats)
print(results)

results_df <- as.data.frame(t(results))  # 转置：每行为变量，每列为分位数
results_df <- round(results_df, 2)  

results_df$PAF_95CI_mort <- paste("(", results_df$Lower95CI, ",", results_df$Upper95CI, ")", sep="")

install.packages("tibble")
library(tibble)

# 将行名变为第一列，列名自定义为 "Variable"（你也可以改为其他）
results_df <- rownames_to_column(results_df, var = "Variable")


library(openxlsx)
write.xlsx(results_df, "~/Desktop/123.xlsx")







#==========================总体PAF的折线图（横坐标是年龄段）==============================================
library(ggplot2)
library(dplyr)

data <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/某一危险因素所有癌症PAF(所有癌种)及case_divided_age.xlsx")

data <- data %>% filter(kind == 0)
data <- data %>% filter(sex == 0)

data$PAF_inc <-data$PAF_inc * 100
data$PAF_lower_inc <-data$PAF_lower_inc * 100
data$PAF_upper_inc <-data$PAF_upper_inc * 100

data <- data %>%
  mutate(kind = case_when(
    kind == 1 ~ "Urban",
    kind == 2 ~ "Rural",
    kind == 0 ~ "Both kind",
    TRUE ~ as.character(kind)
  ))


data1 <- data %>% filter(exposure == "Alcohol"|exposure == "Smoking"|exposure =="Human papillomavirus")
data2 <-data %>% filter(exposure == "Diabetes"|exposure == "Excess bodyweight"|exposure == "Second-hand smoking")
data3 <-data %>% filter(exposure == "Epstein-Barr virus"|exposure == "Hepatitis B virus"|exposure == "Helicobacter pylori")
data4 <-data %>% filter(exposure == "Red meat"|exposure == "Low fruit intake"|exposure == "Low vegetable intake")
data5 <-data %>% filter(exposure == "Clonorchis sinensis"|exposure == "Hepatitis C virus"|exposure == "Physical inactivity"|exposure == "Helicobacter pylori")


ggplot(data5, aes(x = age, y = PAF_inc, group = exposure, color = exposure)) +
  geom_line() +  # 绘制折线图
  geom_point() +  # 添加点以突出显示每个数据点
  labs(x = "Age Group", y = "Population attributable fraction", title = NULL) +
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
    axis.text.x = element_text(size = 9, family = "Times New Roman", color = "black"),
    axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
    axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
    axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
    legend.title = element_text(size = 10, family = "Times New Roman", face = "bold"),  # 设置图例标题样式
    legend.text = element_text(size = 8, family = "Times New Roman"),  # 设置图例文本样式
    axis.title = element_text(size = 10, family = "Times New Roman", face = "bold"),  # 设置轴标题样式
    axis.text = element_text(size = 8, family = "Times New Roman"),  # 设置轴文本样式
    panel.grid = element_blank(),  # 去除背景网格
    strip.text = element_text(size = 5),  # 调整分面标题的字体大小
    legend.position = "right"  # 将图例位置设置为右侧
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # y 轴为百分比格式
  scale_x_discrete(expand = c(0, 0.1)) 


#+  # x 轴为离散值
scale_color_manual(
  name = NULL,
  values = c("Rural" = "#1874CD", "Urban" = "#CD0000", "Both kind" = "#707070"),  # 指定颜色
  guide = guide_legend(nrow = length(unique(data$exposure)))  # 设置图例为单列
)

data1 <- data %>%
  filter(exposure == "Smoking"|exposure == "Alcohol" )

data2 <- data %>%
  filter( exposure == "Physical inactivity" | exposure == "Second-hand smoking")

data3 <- data %>%
  filter(exposure == "Low fruit intake" | exposure == "Low vegetable intake" | exposure == "Red meat" | exposure =="Excess bodyweight"|exposure == "Diabetes"   )

data4 <- data %>%
  filter(exposure == "Hepatitis C virus" | exposure == "Hepatitis B virus" | exposure == "Helicobacter pylori" | exposure == "Epstein-Barr virus")
data5 <- data %>%
  filter(exposure == "Hepatitis C virus" |  exposure == "Helicobacter pylori" | exposure == "Human papillomavirus")


# 查看筛选后的数据框
ggplot(data4, aes(x = age, y = PAF_inc, group = exposure, color = exposure)) +  # 添加 color 映射
  geom_line() +  # 绘制折线图
  geom_point() +  # 添加点以突出显示每个数据点
  labs(x = "Age Group", y = "Population attributable fraction", title = NULL) +
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
    axis.text.x = element_text(size = 9, family = "Times New Roman", color = "black"),
    axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
    axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
    axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
    legend.title = element_text(size = 10, family = "Times New Roman", face = "bold"),  # 设置图例标题样式
    legend.text = element_text(size = 8, family = "Times New Roman"),  # 设置图例文本样式
    axis.title = element_text(size = 10, family = "Times New Roman", face = "bold"),  # 设置轴标题样式
    axis.text = element_text(size = 8, family = "Times New Roman"),  # 设置轴文本样式
    panel.grid = element_blank(),  # 去除背景网格
    strip.text = element_text(size = 5),  # 调整分面标题的字体大小
    legend.position = "right"  # 将图例位置设置为右侧
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # y 轴为百分比格式
  scale_x_discrete(expand = c(0, 0.1)) 


#geom_errorbar(aes(ymin = PAF_lower_inc, ymax = PAF_upper_inc), width = 0.1, color = sex) +  # 添加误差线

ggsave("~/Desktop/广东省肿瘤数据/Plot(2)/总体PAF的折线图/Rural.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 6, 
       height = 3.5, 
       units = "in", 
       dpi = 1200)


#==========================纵坐标site+exposure，横坐标All\sex\kind\age的PAF热图==============================================

install.packages("ggtext")
install.packages("stringr")

library(ggplot2)
library(dplyr)
library(readxl)  
library(ggtext) 
library(stringr)
library(scales)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_age_kind012_sex012.xlsx")
data1<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")
data2<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_all_exposure_ICD10.xlsx")
data3<-read_excel("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_all_exposure_ICD10(all_age).xlsx")




data1 <- data1 %>%
  mutate(age = "All")
data<-rbind(data,data1)
data<-data[,-c(9:14,21:24)]

data3 <- data3 %>%
  mutate(age = "All")
data2<-rbind(data2,data3)
data2<-data2[,-c(14,15)]
data2 <- data2 %>%
  mutate(exposure = "All")

data<-rbind(data,data2)

data$PAF_inc <- as.numeric(data$PAF_inc)
data$PAF_lower_inc <- as.numeric(data$PAF_lower_inc)
data$PAF_upper_inc <- as.numeric(data$PAF_upper_inc)

data$inc_attribute<-data$inc_all*data$PAF_inc
data$inc_attribute_lower<-data$inc_all*data$PAF_lower_inc
data$inc_attribute_upper<-data$inc_all*data$PAF_upper_inc

data$PAF_mort <- as.numeric(data$PAF_mort)
data$PAF_lower_mort <- as.numeric(data$PAF_lower_mort)
data$PAF_upper_mort <- as.numeric(data$PAF_upper_mort)

data$mort_attribute<-data$mort_all*data$PAF_mort
data$mort_attribute_lower<-data$mort_all*data$PAF_lower_mort
data$mort_attribute_upper<-data$mort_all*data$PAF_upper_mort



#data <- data %>%
# mutate(site = str_replace(site, "Gallbladder etc.", "Gallbladder"))

#data$age <- ifelse(data$age == "all", "All", data$age)

data1 <- data %>% filter(kind == 0 & sex ==1 & age == "All")
data2 <- data %>% filter(kind == 0 & sex ==2 & age == "All")
data3 <- data %>% filter(kind == 1 & sex ==0 & age == "All")
data4 <- data %>% filter(kind == 2 & sex ==0 & age == "All")
data5 <- data %>% filter(kind == 0 & sex ==0)
data5 <- data5 %>%
  filter(age == "20-29" | age == "30-39")
data5<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+sex+kind+ICD10+site,data5,sum)
data6 <- data %>% filter(kind == 0 & sex ==0)
data6 <- data6 %>%
  filter(age == "60-69" | age == "70-")
data6<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+sex+kind+ICD10+site,data6,sum)
data7<-data %>% filter(kind == 0 & sex ==0 & age == "All")
data8 <- data %>% filter(kind == 0 & sex ==0)
data8 <- data8 %>%
  filter(age =="40-49" | age == "50-59")
data8<-aggregate(cbind(inc_all,mort_all,inc_attribute,inc_attribute_lower,inc_attribute_upper,mort_attribute,mort_attribute_lower,mort_attribute_upper)~exposure+sex+kind+ICD10+site,data8,sum)


data1<-data1[,-c(2:4,9:14)]
data2<-data2[,-c(2:4,9:14)]
data3<-data3[,-c(2:4,9:14)]
data4<-data4[,-c(2:4,9:14)]
data7<-data7[,-c(2:4,9:14)]
data5<-data5[,-c(2:3)]
data6<-data6[,-c(2:3)]
data8<-data8[,-c(2:3)]
data1 <- data1 %>%
  mutate(factor = "Male")
data2 <- data2 %>%
  mutate(factor = "Female")
data3 <- data3 %>%
  mutate(factor = "Urban")
data4 <- data4 %>%
  mutate(factor = "Rural")
data5 <- data5 %>%
  mutate(factor = "Age=20-39")
data6 <- data6 %>%
  mutate(factor = "Age≥60")
data7 <- data7 %>%
  mutate(factor = "All")
data8 <- data8 %>%
  mutate(factor = "Age=40-59")

data<-rbind(data1,data2,data3,data4,data5,data6,data7,data8)

data$PAF_inc<-data$inc_attribute/data$inc_all
data$PAF_lower_inc<-data$inc_attribute_lower/data$inc_all
data$PAF_upper_inc<-data$inc_attribute_upper/data$inc_all
data$PAF_mort<-data$mort_attribute/data$mort_all
data$PAF_lower_mort<-data$mort_attribute_lower/data$mort_all
data$PAF_upper_mort<-data$mort_attribute_upper/data$mort_all


#填充age20-39的Penis
new_row <- data.frame(
  exposure = "Human papillomavirus",
  site = "Penis",
  sex = 0,
  kind = 0,
  ICD10 = "C60",
  PAF_inc = 0.442,
  PAF_mort = 0.442,
  factor = "Age=20-39",
  stringsAsFactors = FALSE  # 避免将字符向量转换为因子
)
# 添加其他列，这些列的值为 NA
other_columns <- setdiff(names(data), names(new_row))
new_row[other_columns] <- NA

# 将新行添加到原始数据框中
data <- bind_rows(data, new_row)

new_row <- data.frame(
  exposure = "All",
  site = "Penis",
  sex = 0,
  kind = 0,
  ICD10 = "C60",
  PAF_inc = 0.442,
  PAF_mort = 0.442,
  factor = "Age=20-39",
  stringsAsFactors = FALSE  # 避免将字符向量转换为因子
)
# 添加其他列，这些列的值为 NA
other_columns <- setdiff(names(data), names(new_row))
new_row[other_columns] <- NA

# 将新行添加到原始数据框中
data <- bind_rows(data, new_row)


#
data$PAF_inc<-data$PAF_inc*100
data$PAF_mort<-data$PAF_mort*100

data <- data %>%
  mutate(
    site = factor(site, levels = c(
      "Penis",  "Prostate", "Vulva, Vagina", "Breast", "Ovary", "Corpus uteri", "Cervix uteri",
      "Non-hodgkin lymphoma", "Hodgkin disease", "Myeloid leukaemia", "Thyroid", "Bladder", "Kidney",
      "Gallbladder", "Liver", "Pancreas", "Colorectum", "Stomach", "Esophagus", "Traches,bronchus and lung",
      "Larynx", "Nasopharynx", "Oral cavity, pharynx"
    )),  # 设置 site 排序
    
    exposure = factor(exposure, levels = c(
      "Clonorchis sinensis","Hepatitis C virus", "Hepatitis B virus", "Helicobacter pylori", "Human papillomavirus", 
      "Epstein-Barr virus", "Diabetes", "Excess bodyweight", "Red meat", "Low fruit intake", 
      "Low vegetable intake", "Physical inactivity", "Alcohol", "Second-hand smoking", "Smoking","All"
    ))  # 设置 exposure 排序
  )

na_rows <- data[is.na(data$PAF_inc), ]
#data <- na.omit(data)
data$PAF_inc <- ifelse(is.na(data$PAF_inc), 0, data$PAF_inc)
data$PAF_mort <- ifelse(is.na(data$PAF_mort), 0, data$PAF_mort)
colors <- c("#3B5DA0","#3F62A4","#4970AF","#5783B9","#6697C3","#A4C3DD","#CADEEC","#EBF3F8","#FFFFFF","#FDF7DA","#FBE99D","#FAD789","#F8BB6A","#E7995C","#E18D57","#C14336","#AB1F2A")
colors <- c("#567BA9", "white", "#B74A4A")  # 从蓝色到白色再到红色

# PAF_inc
# 转换为均匀分布的值，并将 0 转换为 NA
data <- data %>%
  mutate(PAF_uniform = (rank(PAF_inc) - 1) / (n() - 1),
         PAF_uniform = ifelse(PAF_inc == 0, NA, PAF_uniform))

ggplot(data, aes(x = factor, y = interaction(exposure, site, drop = TRUE), fill = ifelse(is.na(PAF_uniform), NA, PAF_uniform))) +
  geom_tile(color = "white", width = 1, height = 1) +  # 绘制热图格子，格子之间有白色边框
  geom_text(aes(label = ifelse(!is.na(PAF_inc) & PAF_inc != 0, sprintf("%.1f", PAF_inc), "")), family = "Times New Roman", color = "black", size = 2.5) +  # 在每个格子中加入PAF值
  scale_fill_gradientn(
    colors = colors,  # 使用定义的颜色数组
    values = scales::rescale(seq(0, 1, length.out = 100)),  # 使用100个色阶进行过渡
    na.value = "grey",  # 当值为 NA 时填充灰色
    guide = guide_colorbar(
      barwidth = 8,  # 图例色条的宽度
      barheight = 0.4,  # 图例色条的高度
      ticks = FALSE,  # 不显示刻度线
      label.position = "bottom",  # 标签放在色条下方
      labels = percent_format(scale = 1)(rescale(range(data$PAF_inc, na.rm = TRUE))),  # 使用 PAF_inc 的实际值范围
      direction = "horizontal"  # 设置图例为水平排列
    )
  ) +
  labs(x = NULL, y = NULL, fill = "PAF") +  # 添加图例标题
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Times New Roman", color = "black"),
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black"),
    axis.line.x = element_line(color = "grey", linewidth = 0.3),
    axis.line.y = element_line(color = "grey", linewidth = 0.3),
    axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
    axis.ticks = element_line(color = "grey", linewidth = 0.3),
    legend.title = element_text(size = 7, family = "Times New Roman", face = "bold"),  # 设置图例标题样式
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.background = element_rect(color = NA),
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_blank(),  # 隐藏图例文本标签
    panel.grid = element_blank(),  # 去除背景网格
    strip.text = element_text(size = 5),  # 调整分面标题的字体大小
    legend.position = "bottom",  # 将图例位置设置为底部
    panel.background = element_rect(fill = "grey90"),  # 设置坐标轴内背景为灰色
    plot.background = element_rect(fill = "white", color = NA),  # 去除整个图形的边框
    plot.margin = margin(0, 0, 0, 0)  # 移除图形的外边距
  ) +
  scale_x_discrete(limits = c("All", "Male", "Female", "Urban", "Rural", "Age=20-39", "Age=40-59", "Age≥60"), expand = c(0, 0)) +
  scale_y_discrete(
    labels = function(x) {
      sapply(x, function(label) {
        parts <- unlist(strsplit(label, "\\."))
        if (!is.na(parts[1]) && parts[1] == "All") {
          return(paste0("<b>", parts[2], "</b>"))  # 将 site 部分加粗
        } else {
          return(parts[1])
        }
      })
    }
  ) +
  theme(
    axis.text.y = element_markdown()  # 使用 ggtext 来支持 HTML 标签
  )

# PAF_mort
data <- data %>%
  mutate(PAF_uniform = (rank(PAF_mort) - 1) / (n() - 1),
         PAF_uniform = ifelse(PAF_mort == 0, NA, PAF_uniform))

ggplot(data, aes(x = factor, y = interaction(exposure, site, drop = TRUE), fill = ifelse(is.na(PAF_uniform), NA, PAF_uniform))) +
  geom_tile(color = "white", width = 1, height = 1) +  # 绘制热图格子，格子之间有白色边框
  geom_text(aes(label = ifelse(!is.na(PAF_mort) & PAF_inc != 0, sprintf("%.1f", PAF_mort), "")), family = "Times New Roman", color = "black", size = 2.5) +  # 在每个格子中加入PAF值
  scale_fill_gradientn(
    colors = colors,  # 使用定义的颜色数组
    values = rescale(seq(0, 1, length.out = length(colors))),  # 均匀分布的值
    na.value = "grey",  # 当值为 NA 时填充灰色
    guide = guide_colorbar(
      barwidth = 8,  # 图例色条的宽度
      barheight = 0.4,  # 图例色条的高度
      ticks = FALSE,  # 不显示刻度线
      label.position = "bottom",  # 标签放在色条下方
      labels = percent_format(scale = 1)(rescale(range(data$PAF_mort, na.rm = TRUE))),  # 使用 PAF_inc 的实际值范围
      direction = "horizontal"  # 设置图例为水平排列
    )
  ) +
  labs(x = NULL, y = NULL, fill = "PAF") +  # 添加图例标题
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Times New Roman", color = "black"),
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black"),
    axis.line.x = element_line(color = "grey", linewidth = 0.3),
    axis.line.y = element_line(color = "grey", linewidth = 0.3),
    axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
    axis.ticks = element_line(color = "grey", linewidth = 0.3),
    legend.title = element_text(size = 7, family = "Times New Roman", face = "bold"),  # 设置图例标题样式
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.background = element_rect(color = NA),
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_blank(),  # 隐藏图例文本标签
    panel.grid = element_blank(),  # 去除背景网格
    strip.text = element_text(size = 5),  # 调整分面标题的字体大小
    legend.position = "bottom",  # 将图例位置设置为底部
    panel.background = element_rect(fill = "grey90"),  # 设置坐标轴内背景为灰色
    plot.background = element_rect(fill = "white", color = NA),  # 去除整个图形的边框
    plot.margin = margin(0, 0, 0, 0)  # 移除图形的外边距
  ) +
  scale_x_discrete(limits = c("All", "Male", "Female", "Urban", "Rural", "Age=20-39", "Age=40-59", "Age≥60"), expand = c(0, 0)) +
  scale_y_discrete(
    labels = function(x) {
      sapply(x, function(label) {
        parts <- unlist(strsplit(label, "\\."))
        if (!is.na(parts[1]) && parts[1] == "All") {
          return(paste0("<b>", parts[2], "</b>"))  # 将 site 部分加粗
        } else {
          return(parts[1])
        }
      })
    }
  ) +
  theme(
    axis.text.y = element_markdown()  # 使用 ggtext 来支持 HTML 标签
  )

ggsave("~/Desktop/广东省肿瘤数据/Plot(5)/PAF热图-colon和rectum合起来/mort(3age).tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 3.5, 
       height = 14, 
       units = "in", 
       dpi = 1200)




#====================================纵坐标site+exposure，横坐标age的PAF热图=====================================
data <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_age_kind012_sex012.xlsx")
data1<-read_excel("~/Desktop/广东省肿瘤数据/PAF2019/attribute_all_age.xlsx")
data2<-read_excel("~/Desktop/广东省肿瘤数据/PAF2019/attribute_all_exposure_ICD10.xlsx")
data3<-read_excel("~/Desktop/广东省肿瘤数据/PAF2019/attribute_all_exposure_ICD10(all_age).xlsx")


data1 <- data1 %>%
  mutate(age = "All")
data<-rbind(data,data1)
data<-data[,c(1:8,15:20)]

data3 <- data3 %>%
  mutate(age = "All")
data2<-rbind(data2,data3)
data2<-data2[,-c(14,15)]
data2 <- data2 %>%
  mutate(exposure = "All")

data<-rbind(data,data2)

data <- data %>% filter(kind == 1)
data <- data %>% filter(sex == 0)


data <- data %>%
  mutate(
    site = factor(site, levels = c(
      "Penis", "Vulva, Vagina", "Prostate", "Breast", "Ovary", "Corpus uteri", "Cervix uteri",
      "Non-hodgkin lymphoma", "Hodgkin disease", "Myeloid leukaemia", "Thyroid", "Bladder", "Kidney",
      "Gallbladder", "Liver", "Pancreas", "Rectum", "Colon", "Stomach", "Esophagus", "Traches,bronchus and lung",
      "Larynx", "Nasopharynx", "Oral cavity, pharynx"
    )),  # 设置 site 排序
    
    exposure = factor(exposure, levels = c(
      "Hepatitis C virus", "Hepatitis B virus", "Helicobacter pylori", "Human papillomavirus", 
      "Epstein-Barr virus", "Diabetes", "Excess bodyweight", "Red meat", "Low fruit intake", 
      "Low vegetable intake", "Physical inactivity", "Alcohol", "Second-hand smoking", "Smoking","All"
    ))  # 设置 exposure 排序
  )

na_rows <- data[is.na(data$PAF_inc), ]
#data <- na.omit(data)
data$PAF_inc <- ifelse(is.na(data$PAF_inc), 0, data$PAF_inc)
colors <- c("#3B5DA0","#3F62A4","#4970AF","#5783B9","#6697C3","#A4C3DD","#CADEEC","#EBF3F8","#FFFFFF","#FDF7DA","#FBE99D","#FAD789","#F8BB6A","#E7995C","#E18D57","#C14336","#AB1F2A")
colors <- c("#3B5DA0", "white", "#AB1F2A")  # 从蓝色到白色再到红色
# 将 PAF_inc 转换为均匀分布的值，并将 0 转换为 NA
data <- data %>%
  mutate(PAF_uniform = (rank(PAF_inc) - 1) / (n() - 1),
         PAF_uniform = ifelse(PAF_inc == 0, NA, PAF_uniform))

# 绘制热图
ggplot(data, aes(x = age, y = interaction(exposure, site, drop = TRUE), fill = PAF_uniform)) +
  geom_tile(color = "white", width = 1, height = 1) +  # 绘制热图格子，格子之间有白色边框
  geom_text(aes(label = ifelse(!is.na(PAF_inc) & PAF_inc != 0.0, sprintf("%.1f", PAF_inc), "")), family = "Times New Roman", color = "black", size = 2.5) +  # 在每个格子中加入PAF值
  scale_fill_gradientn(
    colors = colors,  # 使用定义的颜色数组
    values = scales::rescale(seq(0, 1, length.out = 100)),  # 使用100个色阶进行过渡
    na.value = "grey",  # 当值为 NA 或 0 时不填充颜色
    guide = guide_colorbar(
      barwidth = 5,  # 图例色条的宽度
      barheight = 0.4,  # 图例色条的高度
      ticks = FALSE,  # 不显示刻度线
      label.position = "bottom",  # 标签放在色条下方
      labels = percent_format(scale = 1)(rescale(range(data$PAF_inc, na.rm = TRUE))),  # 使用 PAF_inc 的实际值范围
      direction = "horizontal"  # 设置图例为水平排列
    )
  ) +
  labs(x = "Age Group", y = NULL, fill = "PAF") +  # 添加图例标题
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5, family = "Times New Roman", color = "black"),
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black"),
    axis.line.x = element_line(color = "grey", linewidth = 0.3),
    axis.line.y = element_line(color = "grey", linewidth = 0.3),
    axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
    axis.ticks = element_line(color = "grey", linewidth = 0.3),
    legend.title = element_text(size = 7, family = "Times New Roman", face = "bold"),  # 设置图例标题样式
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.background = element_rect(color = NA),
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_blank(),  # 隐藏图例文本标签
    panel.grid = element_blank(),  # 去除背景网格
    strip.text = element_text(size = 5),  # 调整分面标题的字体大小
    legend.position = "bottom"  # 将图例位置设置为底部
  ) +
  scale_x_discrete(limits = c("All", "20-29", "30-39", "40-49", "50-59", "60-69", "70-"), expand = c(0, 0)) +
  scale_y_discrete(
    labels = function(x) {
      sapply(x, function(label) {
        parts <- unlist(strsplit(label, "\\."))
        if (!is.na(parts[1]) && parts[1] == "All") {
          return(paste0("<b>", parts[2], "</b>"))  # 将 site 部分加粗
        } else {
          return(parts[1])
        }
      })
    }
  ) +
  theme(
    axis.text.y = element_markdown()  # 使用 ggtext 来支持 HTML 标签
  )


ggsave("~/Desktop/广东省肿瘤数据/Plot(2)/PAF热图/Urban_sex.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 3.5, 
       height = 14, 
       units = "in", 
       dpi = 1200)



rgb_to_hex <- function(r, g, b) {
  sprintf("#%02X%02X%02X", r, g, b)
}

blue_hex <- rgb_to_hex(48, 50, 125)
red_hex <- rgb_to_hex(142, 0, 40)
color_1 <- rgb_to_hex(189, 60, 51)
color_2 <- rgb_to_hex(228, 148, 90)
color_3 <- rgb_to_hex(249, 198, 118)
color_4 <- rgb_to_hex(251, 236, 171)
color_5 <- rgb_to_hex(255, 255, 255)
color_6 <- rgb_to_hex(209, 226, 239)
color_7 <- rgb_to_hex(209, 226, 239)
color_8 <- rgb_to_hex(127, 169, 205)
color_9 <- rgb_to_hex(82, 124, 180)

#colors <- c(blue_hex, color_9, color_8, color_7, color_6, color_5, color_4, color_3, color_2, color_1, red_hex)

# 归一化每一行的PAF
data <- data %>%
  group_by(exposure, site) %>%
  mutate(PAF_rank = rank(-PAF_inc, ties.method = "first")) %>%  # 计算每行的PAF排名
  ungroup()  # 取消分组，以便后续操作不受分组影响

ggplot(data, aes(x = age, y = interaction(exposure, site, drop = TRUE), fill = PAF_rank)) +
  geom_tile(color = "white", width = 1, height = 1) +  # 绘制热图格子，格子之间有白色边框
  geom_text(aes(label = ifelse(!is.na(PAF_inc) & PAF_inc != 0.0, sprintf("%.1f", PAF_inc), "")), family = "Times New Roman", color = "black", size = 2.5) +  # 在每个格子中加入PAF值
  scale_fill_gradientn(
    colors = c("seagreen", "yellow", "red2"),  # 自定义渐变色条
    values = scales::rescale(c(7, 4, 1)),  # 设置排名范围对应的渐变色条
    guide = guide_colorbar(
      barwidth = 5,  # 图例色条的宽度
      barheight = 0.4,  # 图例色条的高度
      ticks = FALSE,  # 不显示刻度线
      label.position = "bottom",  # 标签放在色条下方
      labels = c("2", "6"),  # 只显示1和6的标签
      direction = "horizontal"  # 设置图例为水平排列
    )
  ) +
  labs(x = "Age Group", y = NULL, fill = "Rank of PAF") +  # 添加图例标题
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5, family = "Times New Roman", color = "black"),
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black"),
    axis.line.x = element_line(color = "grey", linewidth = 0.3),
    axis.line.y = element_line(color = "grey", linewidth = 0.3),
    axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
    axis.ticks = element_line(color = "grey", linewidth = 0.3),
    legend.title = element_text(size = 7, family = "Times New Roman", face = "bold"),  # 设置图例标题样式
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.background = element_rect(color = NA),
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_text(size = 6, family = "Times New Roman", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),  # 去除背景网格
    strip.text = element_text(size = 5),  # 调整分面标题的字体大小
    plot.margin = margin(0, 0, 0, 0),  # 去掉周围的额外边距
    legend.position = "bottom"  # 将图例位置设置为底部
  ) +
  scale_x_discrete(limits = c("All", "20-29", "30-39", "40-49", "50-59", "60-69", "70-"), expand = c(0, 0)) +
  scale_y_discrete(
    labels = function(x) {
      sapply(x, function(label) {
        parts <- unlist(strsplit(label, "\\."))
        if (!is.na(parts[1]) && parts[1] == "All") {
          return(paste0("<b>", parts[2], "</b>"))  # 将 site 部分加粗
        } else {
          return(parts[1])
        }
      })
    }
  ) +
  theme(
    axis.text.y = element_markdown()  # 使用 ggtext 来支持 HTML 标签
  )



#==========================按癌种区分的各种危险因素PAF百分比图==============================================

install.packages("readxl")
install.packages("dplyr")
install.packages("ggplot2")

library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)

# 读取数据
df <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_all_age.xlsx")

df <- df %>% filter(kind == 0)
df <- df %>% filter(sex == 0)

df_complete_na <- df %>% filter(is.na(site))
df_complete_na_all <- df[!complete.cases(df), ]

# 计算每个癌种的总病例数
df_sum <- df %>%
  group_by(site) %>%
  mutate(total_cases = sum(inc_attribute)) %>%
  ungroup() %>%
  # 计算每个暴露因素占总病例数的百分比
  mutate(percentage = (inc_attribute / total_cases) * 100)

# 定义暴露因素的顺序
exposure_order <- c(
  "Smoking", "Second-hand smoking", "Alcohol", "Physical inactivity",  # 行为因子的顺序
  "Low fruit intake", "Low vegetable intake", "Red meat",  # 饮食因子的顺序
  "Diabetes", "Excess bodyweight" , # 代谢因子的顺序
  "Human papillomavirus", "Epstein-Barr virus", "Helicobacter pylori", "Hepatitis B virus" , "Hepatitis C virus"  # 感染因子的顺序
  
)

# 确保每个癌种都包含所有暴露因素，并按暴露因素的顺序排列
df_complete <- df_sum %>%
  complete(site, exposure = exposure_order, fill = list(inc_attribute = 0, percentage = 0)) %>%
  mutate(
    exposure = factor(exposure, levels = exposure_order)  # 确保暴露因子按顺序排列
  )

# 将暴露因素的颜色映射到指定的颜色
factor_colors <- c(
  "Smoking" = "#AB1F2A", 
  "Second-hand smoking" ="#C14336" , 
  "Alcohol" = "#E18D57", 
  "Physical inactivity" = "#E7995C",
  "Low fruit intake" = "#F8BB6A", 
  "Low vegetable intake" = "#FAD789", 
  "Red meat" = "#FBE99D",
  "Diabetes" = "#FDF7DA", 
  "Excess bodyweight" = "#A4C3DD",
  "Human papillomavirus" = "#6697C3", 
  "Epstein-Barr virus" = "#5783B9", 
  "Helicobacter pylori" = "#4970AF", 
  "Hepatitis B virus" = "#3F62A4", 
  "Hepatitis C virus"= "#3B5DA0"
)

# 定义癌种的顺序
site_order <- c(
  "Penis", "Vulva, Vagina", "Prostate", "Breast", "Ovary", "Corpus uteri", "Cervix uteri",  "Non-hodgkin lymphoma", "Hodgkin disease","Myeloid leukaemia", "Thyroid","Bladder", "Kidney", "Gallbladder", "Liver", "Pancreas", "Rectum", "Colon", "Stomach","Esophagus", "Traches,bronchus and lung", "Larynx", "Nasopharynx", "Oral cavity, pharynx"
)

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
df_complete <- df_complete %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序



labels <- c(
  "Helicobacter pylori" = expression(italic("Helicobacter pylori")), 
  "Smoking" = "Smoking", 
  "Second-hand smoking" = "Second-hand smoking", 
  "Alcohol" = "Alcohol", 
  "Physical inactivity" = "Physical inactivity",
  "Low fruit intake" = "Low fruit intake", 
  "Low vegetable intake" = "Low vegetable intake", 
  "Red meat" = "Red meat",
  "Human papillomavirus" = "Human papillomavirus", 
  "Epstein-Barr virus" = "Epstein-Barr virus", 
  "Hepatitis B virus" = "Hepatitis B virus", 
  "Hepatitis C virus" = "Hepatitis C virus",
  "Diabetes" = "Diabetes", 
  "Excess bodyweight" = "Excess bodyweight"
)

ggplot(df_complete, aes(x = percentage, y = site, fill = exposure)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = factor_colors, labels = labels) +  # 设置颜色和标签
  labs(
    x = "Proportion of cancer cases attributable to modifiable risk factors",  # 横坐标标题
    y = NULL  # 去掉纵坐标标题
  ) +
  theme(
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
    axis.text.x = element_text(size = 9, family = "Times New Roman", color = "black"),
    axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
    axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
    axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
    legend.position = "right",  # 图例放置在右侧
    legend.key.size = unit(0.4, "cm"),  # 缩小图例的大小
    legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),
    panel.grid = element_blank(),  # 去掉网格线
    panel.background = element_blank(),  # 去掉背景色
    plot.background = element_blank(),  # 去掉整个图表背景色
    axis.title.x = element_text(size = 10, family = "Times New Roman", face = "bold", color = "black")  # 横坐标标题字体大小
  ) +
  scale_x_continuous(labels = scales::percent_format(scale = 1), expand = c(0, 0)) +  # x 轴为百分比格式
  guides(fill = guide_legend(title = NULL))  # 去掉图例标题

ggsave("~/Desktop/广东省肿瘤数据/Plot(2)/按癌种区分的各种危险因素PAF百分比图/Both_kind_sex.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 8, 
       height = 3.5, 
       units = "in", 
       dpi = 1000)


#============================按癌种区分的所有危险因素PAF条图==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)


data <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_case_all_exposure_ICD10(all_age).xlsx")

data1 <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/所有危险因素-所有癌症PAF及case.xlsx")
data1 <- data1 %>%
  mutate(site = "All")
data1 <- data1 %>%
  mutate(ICD10 = "All")

data<-rbind(data,data1)

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 0)

data$PAF_inc<-data$PAF_inc*100
data$PAF_lower_inc<-data$PAF_lower_inc*100
data$PAF_upper_inc<-data$PAF_upper_inc*100
data$PAF_inc1 <- sprintf("%0.2f", data$PAF_inc)

data$PAF_mort<-data$PAF_mort*100
data$PAF_lower_mort<-data$PAF_lower_mort*100
data$PAF_upper_mort<-data$PAF_upper_mort*100
data$PAF_mort1 <- sprintf("%0.2f", data$PAF_mort)

data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))


# 将癌种按定义的顺序排序，并设置为因子
combined_data <- data %>%
  mutate(
    sex = factor(sex, levels = c("Total","Male", "Female"))
  )
combined_data <- combined_data %>%
  group_by(sex) %>%
  arrange(PAF_inc, .by_group = TRUE) %>%  # 对每个 sex 组内按 PAF_mort 排序
  mutate(site = factor(site, levels = unique(site))) %>%  # 确保排序后不改变 site 的因子级别
  ungroup() 
#facet_grid(rows = vars(sex),  cols = NULL, scales = "free_y", , space = "free_y", switch = "y", shrink = TRUE, facets = NULL,margins = FALSE, drop = TRUE, as.table = TRUE)
ggplot(combined_data, aes(x = PAF_inc, y = site, fill = sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +  # 使用 dodge 以确保条形图不重叠
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100),  labels = scales::percent_format(scale = 1), expand = c(0, 0)) +  # 扩展x轴范围
  scale_y_discrete() +  # 使用 discrete scale for y axis
  scale_fill_manual(values = c("Male" = "#3F62A4", "Female" = "#AB1F2A", "Total" = "#707070")) +  # 男性为蓝色，女性为红色，总计为绿色
  labs(x =  NULL, y = NULL) +
  facet_grid(rows = vars(sex), scales = "free_y", space = "free_y", switch = "y") +  # 将面板标签移动到左侧
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 9, family = "Times New Roman", color = "black"),
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black"),
    axis.line.x = element_line(color = "grey", size = 0.3),
    axis.line.y = element_line(color = "grey", size = 0.3),
    axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
    axis.ticks = element_line(color = "grey", size = 0.3),
    legend.position = "bottom",  # 将图例放置在右侧
    legend.title = element_blank(),  # 删除图例标题
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_text(size = 9, family = "Times New Roman", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
    strip.placement = "left"  # 将facet标签放置在左侧
  ) +
  # 添加inc_attribute文本标签到x = 115的位置
  geom_text(aes(x = 115 , y = site, label = inc_95CI), 
            hjust = 0, size = 2.5, family = "Times New Roman", color = "black") +  # 将标签放到x = 115位置
  coord_cartesian(clip = "off") + 
  # 只在'Male'面板上添加"Attributable Cases"标签
  geom_text(data = subset(combined_data, sex == "Total"), 
            aes(x = 150, y = 25.5, label = "Attributable Cases"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Times New Roman", color = "black")+ 
  # 只在'Male'面板上添加"PAF(%)"标签
  geom_text(data = subset(combined_data, sex == "Total"), 
            aes(x = 60, y = 25.5, label = "PAF(%)"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Times New Roman", color = "black")+
  # 在条形图右侧添加总和标签
  geom_text(data = data, aes(x = PAF_inc+1, y = site, label = PAF_inc1), 
            vjust = 0.5, hjust = 0, size = 2.5, family = "Times New Roman", color = "black")

#geom_text(data = subset(combined_data, sex == "Male"), 
#aes(x = 60, y = 21, label = "PAF(%)"),
#hjust = 1.1, vjust = 1.1, size = 3.5, family = "Times New Roman", fontface = "bold", color = "black")
#geom_errorbar(aes(xmin = PAF_lower_inc, xmax = PAF_upper_inc), width = 0.3, color = "black", position = position_dodge(0.8)) +

#==================================最新：按癌种区分的所有危险因素PAF条图====================================
library(gridExtra)
library(ggplot2)
library(grid)
library(readxl)
library(dplyr)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/attribute_case_all_exposure_ICD10(all_age).xlsx")

data1 <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/所有危险因素-所有癌症PAF及case.xlsx")
data1 <- data1 %>%
  mutate(site = "All")
data1 <- data1 %>%
  mutate(ICD10 = "All")

data<-rbind(data,data1)

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 1)

data$PAF_inc<-data$PAF_inc*100
data$PAF_lower_inc<-data$PAF_lower_inc*100
data$PAF_upper_inc<-data$PAF_upper_inc*100
data$PAF_inc1 <- sprintf("%0.1f", data$PAF_inc)

data$PAF_mort<-data$PAF_mort*100
data$PAF_lower_mort<-data$PAF_lower_mort*100
data$PAF_upper_mort<-data$PAF_upper_mort*100
data$PAF_mort1 <- sprintf("%0.1f", data$PAF_mort)

data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))
#"#F9AE78","#ffd36a","#ffc99a","#e2d5ba","#5176A6","#F26B8F","#ffd899","#6abac4","#fa8072","#478ecc""#D44E6B""#b395bd""#fbd178""#abd9d0","#74a9c5", "Female" = "#dd7389" , "Total" = "#f6cd96"
data1 <- data %>% filter(sex == "Total")
data1 <- data1 %>%
  arrange(PAF_inc) %>%  # 对每个 sex 组内按 PAF_mort 排序
  mutate(site = factor(site, levels = unique(site))) # 确保排序后不改变 site 的因子级别
#"Male" = "#74a9c5""#3366CC", "Female" = "#dd7389""#FF99CC" , "Total" = "#D0D4D7""#f6cd96"
plot1 <- ggplot(data1, aes(x = PAF_inc, y = site)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, fill = "#D0D4D7") +  # 设置条形颜色为绿色
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100), labels = scales::percent_format(scale = 1), expand = c(0, 0)) +  # 扩展x轴范围
  scale_y_discrete() +  # 使用 discrete scale for y axis
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),  # 去除x轴标签
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black"),  # 纵坐标标签样式
    axis.line.x = element_line(color = "grey", size = 0.3),
    axis.line.y = element_line(color = "grey", size = 0.3),
    axis.ticks = element_line(color = "grey", size = 0.3),
    legend.position = "bottom",  # 将图例放置在底部
    legend.title = element_blank(),  # 删除图例标题
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_text(size = 9, family = "Times New Roman", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
    strip.placement = "left",  # 将facet标签放置在左侧
    plot.margin = margin(l = 50)  # 调整左边距以腾出空间给文本标签
  ) +
  # 添加inc_attribute文本标签到x = 115的位置
  geom_text(aes(x = 115 , y = site, label = inc_95CI), 
            hjust = 0, size = 2.5, family = "Times New Roman", color = "black") +  # 将标签放到x = 115位置
  coord_cartesian(clip = "off") +
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = PAF_inc + 1, y = site, label = PAF_inc1), 
            vjust = 0.5, hjust = 0, size = 2.5, family = "Times New Roman", color = "black") +
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Total", gp = gpar(fontsize = 10, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -107, xmax = 0, ymin = 0, ymax = length(unique(data1$site))  # 设置文本的坐标
  )+ 
  # 只在'Male'面板上添加"Attributable Cases"标签
  geom_text(data = data1, 
            aes(x = 147, y = 25, label = "Attributable Cases"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Times New Roman", color = "black")+ 
  # 只在'Male'面板上添加"PAF(%)"标签
  geom_text(data = data1, 
            aes(x = 60, y = 25, label = "PAF(%)"),
            hjust = 1.1, vjust = 1.1, size = 3.5, family = "Times New Roman", color = "black")


data2 <- data %>% filter(sex == "Male")
data2 <- data2 %>%
  arrange(PAF_inc) %>%  # 对每个 sex 组内按 PAF_mort 排序
  mutate(site = factor(site, levels = unique(site))) # 确保排序后不改变 site 的因子级别

plot2 <- ggplot(data2, aes(x = PAF_inc, y = site)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, fill = "#3366CC") +  # 设置条形颜色为绿色
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100), labels = scales::percent_format(scale = 1), expand = c(0, 0)) +  # 扩展x轴范围
  scale_y_discrete() +  # 使用 discrete scale for y axis
  labs(x = NULL, y = NULL) +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),  # 去除x轴标签
    axis.title.x = element_blank(),  # 去除x轴标题
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black"),  # 纵坐标标签样式
    axis.line.x = element_line(color = "grey", size = 0.3),
    axis.line.y = element_line(color = "grey", size = 0.3),
    axis.ticks = element_line(color = "grey", size = 0.3),
    legend.position = "bottom",  # 将图例放置在底部
    legend.title = element_blank(),  # 删除图例标题
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_text(size = 9, family = "Times New Roman", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
    strip.placement = "left",  # 将facet标签放置在左侧
    plot.margin = margin(l = 50)  # 调整左边距以腾出空间给文本标签
  ) +
  # 添加inc_attribute文本标签到x = 115的位置
  geom_text(aes(x = 115 , y = site, label = inc_95CI), 
            hjust = 0, size = 2.5, family = "Times New Roman", color = "black") +  # 将标签放到x = 115位置
  coord_cartesian(clip = "off") +
  # 在条形图右侧添加总和标签
  geom_text(data = data2, aes(x = PAF_inc + 1, y = site, label = PAF_inc1), 
            vjust = 0.5, hjust = 0, size = 2.5, family = "Times New Roman", color = "black") +
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Male", gp = gpar(fontsize = 10, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -107, xmax = 0, ymin = 0, ymax = length(unique(data1$site))  # 设置文本的坐标
  )

data3 <- data %>% filter(sex == "Female")
data3 <- data3 %>%
  arrange(PAF_inc) %>%  # 对每个 sex 组内按 PAF_mort 排序
  mutate(site = factor(site, levels = unique(site))) # 确保排序后不改变 site 的因子级别

plot3 <- ggplot(data3, aes(x = PAF_inc, y = site)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, fill = "#FF99CC" ) +  # 设置条形颜色为灰色
  scale_x_continuous(
    limits = c(0, 150), 
    breaks = c(0, 50, 100), 
    labels = c("0", "50", "100"),  # 设置横坐标轴刻度标签
    expand = c(0, 0)
  ) +  # 扩展x轴范围并设置刻度标签
  scale_y_discrete() +  # 使用 discrete scale for y axis
  labs(x = NULL, y = NULL) +  # 删除x轴标题
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 9, family = "Times New Roman", color = "black"),  # 横坐标标签样式
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black"),  # 纵坐标标签样式
    axis.line.x = element_line(color = "grey", size = 0.3),
    axis.line.y = element_line(color = "grey", size = 0.3),
    axis.ticks = element_line(color = "grey", size = 0.3),
    legend.position = "bottom",  # 将图例放置在底部
    legend.title = element_blank(),  # 删除图例标题
    legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
    legend.margin = margin(t = -5),  # 设置负的上边距，使图例更接近坐标轴
    legend.text = element_text(size = 9, family = "Times New Roman", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
    strip.placement = "left",  # 将facet标签放置在左侧
    plot.margin = margin(l = 50)  # 调整左边距以腾出空间给文本标签
  ) +
  # 添加inc_attribute文本标签到x = 115的位置
  geom_text(aes(x = 115 , y = site, label = inc_95CI), 
            hjust = 0, size = 2.5, family = "Times New Roman", color = "black") +  # 将标签放到x = 115位置
  coord_cartesian(clip = "off") +
  # 在条形图右侧添加总和标签
  geom_text(data = data3, aes(x = PAF_inc + 1, y = site, label = PAF_inc1), 
            vjust = 0.5, hjust = 0, size = 2.5, family = "Times New Roman", color = "black") +
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 10, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -107, xmax = 0, ymin = 0, ymax = length(unique(data1$site))  # 设置文本的坐标
  )


plot1 <- plot1 + theme(plot.margin = margin(t = 10, r = 10, b = 1, l = 50))  # Adjust for plot1
plot2 <- plot2 + theme(plot.margin = margin(t = 5, r = 10, b = 1, l = 50))  # Adjust for plot2
plot3 <- plot3 + theme(plot.margin = margin(t = 5, r = 10, b = 10, l = 50))  # Increase bottom margin for plot3

# 使用 grid.arrange() 将三个图排列
plot <- grid.arrange(plot1, plot2, plot3, ncol = 1, heights = c(1.12,0.9,1.1))



ggsave("~/Desktop/数据/广东省肿瘤归因/分析结果/Plot/Plot(6.9)/按癌种区分的PAF条图/inc-kind1.tiff", 
       plot = plot, 
       device = "tiff", 
       width = 7.1, 
       height = 9, 
       units = "in", 
       dpi = 1000)



#====================================================================
library(readxl) 
library(ggplot2)
library(dplyr)


data <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_case_all_exposure_ICD10(all_age).xlsx")

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 0)

data$PAF_inc<-data$PAF_inc*100
data$PAF_lower_inc<-data$PAF_lower_inc*100
data$PAF_upper_inc<-data$PAF_upper_inc*100
data$PAF_inc1 <- sprintf("%0.2f", data$PAF_inc)

data$PAF_mort<-data$PAF_mort*100
data$PAF_lower_mort<-data$PAF_lower_mort*100
data$PAF_upper_mort<-data$PAF_upper_mort*100
data$PAF_mort1 <- sprintf("%0.2f", data$PAF_mort)

data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))


# 将癌种按定义的顺序排序，并设置为因子
combined_data <- data %>%
  mutate(
    sex = factor(sex, levels = c("Male", "Female", "Total"))
  )

combined_data <- combined_data %>%
  group_by(sex) %>%
  arrange(PAF_inc, .by_group = TRUE) %>%  # 对每个 sex 组内按 PAF_mort 排序
  mutate(site = factor(site, levels = unique(site))) %>%  # 确保排序后不改变 site 的因子级别
  ungroup()

# 使用 base R 的 barplot 函数

# 为每个面板分别绘制条形图
par(mfrow = c(3, 1), mar = c(2, 10, 1, 5), oma = c(0, 0, 0, 0))  # 设定面板布局（3行1列），缩窄面板之间的边距

# Male面板
male_data <- combined_data_sorted %>% filter(sex == "Male")
barplot(male_data$PAF_inc, names.arg = male_data$site, 
        col = "#3F62A4", las = 2,  # 设置颜色为蓝色，并旋转标签
        cex.names = 0.8, # 缩小横坐标标签
        main = "",  # 不显示标题
        border = NA, 
        horiz = TRUE,  # 水平条形图
        axes = TRUE,  # 显示坐标轴
        xaxt = "n")   # 不显示横坐标标签

# Female面板
female_data <- combined_data_sorted %>% filter(sex == "Female")
barplot(female_data$PAF_inc, names.arg = female_data$site, 
        col = "#AB1F2A", las = 2, 
        cex.names = 0.8, 
        main = "", 
        border = NA, 
        horiz = TRUE,
        axes = TRUE,  # 显示坐标轴
        xaxt = "n")   # 不显示横坐标标签

# Total面板
total_data <- combined_data_sorted %>% filter(sex == "Total")
barplot(total_data$PAF_inc, names.arg = total_data$site, 
        col = c("#3F62A4", "#AB1F2A", "#707070")[factor(total_data$sex)], # 为每个性别设置颜色
        las = 2, 
        cex.names = 0.8, 
        main = "", 
        border = NA, 
        horiz = TRUE, 
        axes = TRUE)  # 显示坐标轴


# 添加图例，放置在横坐标标签下方
legend("right", 
       legend = c("Male", "Female", "Total"), 
       fill = c("#3F62A4", "#AB1F2A", "#707070"), 
       title = "Sex", 
       cex = 0.8, 
       horiz = TRUE, 
       xpd = TRUE)

# 恢复默认设置
par(mfrow = c(1, 1))




#====================================================================
#+
# 使用geom_segment添加一条横跨两个面板的参考线
#geom_segment(
#aes(x = -Inf, xend = Inf, y = 0, yend = 0),  # 横跨整个x轴并在y = 0处画线
#color = "black", size = 0.2, linetype = "solid"
# )


#============================按exposure区分的各癌种case条图==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)


data <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_all_age.xlsx")

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 0)
#



# 定义癌种的顺序
site_order <- c(
  "Oral cavity, pharynx","Nasopharynx","Larynx","Traches,bronchus and lung","Esophagus","Stomach","Colon", "Rectum","Pancreas","Liver","Gallbladder", "Kidney","Bladder","Thyroid", "Myeloid leukaemia","Hodgkin disease","Non-hodgkin lymphoma","Cervix uteri","Corpus uteri","Ovary","Breast","Prostate","Vulva, Vagina","Penis"
)

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
data <- data %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序


data <- data %>%
  group_by(exposure, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(inc_attribute))  #%>%
#ungroup() %>%
#mutate(exposure = factor(exposure, levels = unique(exposure[order(total_cases)])))  # 按total_cases从高到低重新排序


#labels <- c(
"Helicobacter pylori" = expression(italic("Helicobacter pylori")), 
"Smoking" = "Smoking", 
"Second-hand smoking" = "Second-hand smoking", 
"Alcohol" = "Alcohol", 
"Physical inactivity" = "Physical inactivity",
"Low fruit intake" = "Low fruit intake", 
"Low vegetable intake" = "Low vegetable intake", 
"Red meat" = "Red meat",
"Human papillomavirus" = "Human papillomavirus", 
"Epstein-Barr virus" = "Epstein-Barr virus", 
"Hepatitis B virus" = "Hepatitis B virus", 
"Hepatitis C virus" = "Hepatitis C virus",
"Diabetes" = "Diabetes", 
"Excess bodyweight" = "Excess bodyweight"
)


#custom_colors <- c(
#"#AB1F2A", "#B73330", "#C54C3A" ,"#D7764C" ,"#E29058", "#E5965B", "#EDA661", "#F7B969", "#F9C97A", "#FAD88A",
#"#FAE296", "#FBECAA", "#FCF3CC", "#DEE4DB", "#ABC7DC" ,"#86ADD0", "#6596C2", "#5C8ABC", "#547FB7", "#4C74B1",
#"#456BAB" ,"#4063A5", "#3D5FA2", "#3B5DA0")

custom_colors <- c("#AB1F2A", "#B5302F", "#C14537", "#CF6544", "#DC8251", "#EAA363", "#F9CC7D", "#FAD88B", "#FAE195", "#FAEAA6",
                   "#FBF0C1", "#EBEAD4" ,"#C7D7DB", "#9FBFD8", "#80A9CD", "#6495C1", "#5C8BBC", "#5581B8", "#4E77B3", "#486FAD",
                   "#4367A8" ,"#3F61A4" ,"#3C5EA1" ,"#3F62A4")

panel_labels_sex <- c("2" = "Female", 
                      "1" = "Male", 
                      "0" = "Total")

data <- data %>%
  mutate(sex = as.character(sex))

# 定义暴露因素的顺序
exposure_order <- list(
  "Total" = c("Hepatitis C virus","Physical inactivity","Low vegetable intake","Helicobacter pylori",
              "Second-hand smoking","Red meat","Excess bodyweight" ,"Epstein-Barr virus",
              "Human papillomavirus","Hepatitis B virus" ,"Low fruit intake","Diabetes","Alcohol","Smoking"
  ),
  "Male" = c("Hepatitis C virus","Physical inactivity","Low vegetable intake","Helicobacter pylori",
             "Second-hand smoking","Red meat","Excess bodyweight" ,"Epstein-Barr virus",
             "Human papillomavirus","Hepatitis B virus" ,"Low fruit intake","Diabetes","Alcohol","Smoking" 
  ),
  "Female" = c("Hepatitis C virus","Physical inactivity","Low vegetable intake","Helicobacter pylori",
               "Second-hand smoking","Red meat","Excess bodyweight" ,"Epstein-Barr virus",
               "Hepatitis B virus" ,"Low fruit intake","Diabetes","Alcohol","Smoking","Human papillomavirus")
)

# 确保每个 sex 面板中的 exposure 按顺序排列
data <- data %>%
  mutate(exposure = case_when(
    sex == "Total" ~ factor(exposure, levels = exposure_order[["Total"]]),
    sex == "Male" ~ factor(exposure, levels = exposure_order[["Male"]]),
    sex == "Female" ~ factor(exposure, levels = exposure_order[["Female"]]),
    TRUE ~ exposure  # 对于其他 sex 情况不做改变
  ))

data <- data %>% filter(sex == 0)

exposure_order <- c("Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Low vegetable intake","Helicobacter pylori",
                    "Second-hand smoking","Red meat","Excess bodyweight" ,"Epstein-Barr virus",
                    "Human papillomavirus","Hepatitis B virus" ,"Low fruit intake","Diabetes","Alcohol","Smoking"
)

data <- data %>%
  mutate(exposure = factor(exposure, levels = exposure_order))

ggplot(data, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 40000), breaks = c(0, 10000, 20000,30000,40000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = "Exposure", fill = "Cancer Site") + 
  labs(y = NULL) +
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 9, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines")) +
  facet_wrap(~ sex, ncol = 1, scales = "free_y", 
             labeller = as_labeller(panel_labels_sex), 
             strip.position = "left") + # 将面板标签放置在左侧
  # 在条形图右侧添加总和标签
  geom_text(data = data, aes(x = total_cases+500, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")


ggsave("~/Desktop/广东省肿瘤数据/Plot(3)/按exposure区分的各癌种case条图/Both_kind.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 7, 
       height = 3.5, 
       units = "in", 
       dpi = 1000)

#改变横坐标轴大小（城乡合并）
ggplot(data, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 40000), breaks = c(0, 10000,20000,30000, 40000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = "Exposure", fill = "Cancer Site") + 
  labs(y = NULL) +
  scale_fill_viridis(discrete = TRUE, guide = guide_legend(ncol = 1), option = "D") +  # 设置图例为一列
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 9, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank())+
  # 在条形图右侧添加总和标签
  geom_text(data = data, aes(x = total_cases+500, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")

#============================按exposure区分的各癌种case条图inc（更新）==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)
library(grid)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 0)

# 定义癌种的顺序
site_order <- c(
  "Traches,bronchus and lung","Nasopharynx","Oral cavity, pharynx","Larynx",
  "Liver","Colorectum","Stomach","Esophagus","Pancreas","Gallbladder",
  "Cervix uteri","Breast","Corpus uteri","Prostate","Ovary","Vulva, Vagina", "Penis",
  "Thyroid","Bladder","Kidney","Myeloid leukaemia","Hodgkin disease"
)

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
data <- data %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序

# 定义初始的颜色节点
nodes <- c( "#bc6356","#dc917b","#eabaa1","#fee3ce")

# 创建插值函数
color_interpolation <- colorRampPalette(nodes)

# 生成22种颜色
colors1 <- color_interpolation(7)
print(colors1)

#custom_colors <- c("#3581b7","#4A92C3","#5FA3CF","#75B5DC",
"#9B3F5C" ,"#AD5873", "#BF718B", "#D08BA7", "#DFA9C5", "#EBCCE2",
"#BC6356" ,"#CC7968", "#DC917B", "#E3A58E", "#EABAA1", "#F3CEB7" ,"#FEE3CE",
"#4d7e54","#669877","#81b095","#a4cbb7","#cfeadf"
)
custom_colors <- c("#4FA8C7", "#6BB8D5", "#83C4E0", "#9ED8E7", 
                   "#D17D8E", "#D99CA8", "#E0A3B4", "#E7B8C0", "#F1C5D4", "#F5D0E1", 
                   "#D64F3A", "#E7725C", "#EB8A70", "#ED9D85", "#F1B1A3", "#F5C6B9","#F7D0C5",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8"
)
custom_colors <- c("#567CD0","#82A3E0","#A9C2EA" ,"#D6E8FA", 
                   "#F2A7C0", "#F5B6CC", "#F8C5D8", "#FBD4E4", "#FEE3ED", "#FFF2F7", 
                   "#E47A6D", "#E68E80", "#E8A292", "#EAB6A4", "#EDCAB6", "#EFD0A8", "#F1D6B9",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8"
)


data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

data <- data %>%
  mutate(sex = as.character(sex))

data <- data %>%
  group_by(exposure, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(inc_attribute))  #%>%


plot0 <- ggplot(data, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 26000), breaks = c(0,5000, 10000, 15000, 20000, 25000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),  # 设置图例标题
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 设置图例文本样式
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -12000, xmax = -10000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示



data1 <- data %>% filter(sex == "Total")
exposure_order1 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Low vegetable intake","Helicobacter pylori","Second-hand smoking",
  "Red meat","Excess bodyweight","Human papillomavirus","Epstein-Barr virus",
  "Hepatitis B virus", "Low fruit intake","Diabetes", "Alcohol","Smoking"
)
data1 <- data1 %>%
  mutate(exposure = factor(exposure, levels = exposure_order1))  # 按照提供的顺序调整 site 的因子顺序

plot1 <- ggplot(data1, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 26000), breaks = c(0,5000, 10000, 15000, 20000, 25000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = NULL, y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 去除x轴标签
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Total", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -6000, xmax = -4000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示


data2 <- data %>% filter(sex == "Male")
exposure_order2 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Low vegetable intake","Second-hand smoking",
  "Human papillomavirus","Helicobacter pylori","Excess bodyweight","Red meat","Epstein-Barr virus",
  "Diabetes", "Low fruit intake","Hepatitis B virus", "Alcohol","Smoking"
)
data2 <- data2 %>%
  mutate(exposure = factor(exposure, levels = exposure_order2))  # 按照提供的顺序调整 site 的因子顺序

plot2 <- ggplot(data2, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 26000), breaks = c(0,5000, 10000, 15000, 20000, 25000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = "Exposure", fill = "Cancer Site") + 
  labs(y = NULL) +
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_blank(),  # 去除x轴标签
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Male", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -6000, xmax = -4000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data2, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

data3 <- data %>% filter(sex == "Female")
exposure_order3 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Low vegetable intake","Physical inactivity","Helicobacter pylori",
  "Smoking","Red meat","Second-hand smoking","Hepatitis B virus","Alcohol","Epstein-Barr virus",
  "Excess bodyweight","Low fruit intake","Diabetes","Human papillomavirus"
)
data3 <- data3 %>%
  mutate(exposure = factor(exposure, levels = exposure_order3))  # 按照提供的顺序调整 site 的因子顺序

plot3 <- ggplot(data3, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 26000), breaks = c(0,5000, 10000, 15000, 20000, 25000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.text = element_blank(),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -6000, xmax = -4000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data3, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

plot1 <- plot1 + theme(plot.margin = margin(t = 10, r = 100, b = 1, l = 50))  # Adjust for plot1
plot2 <- plot2 + theme(plot.margin = margin(t = 5, r = 100, b = 1, l = 50))  # Adjust for plot2
plot3 <- plot3 + theme(plot.margin = margin(t = 5, r = 100, b = 10, l = 50))  # Increase bottom margin for plot3

# 使用 grid.arrange() 将三个图排列
plot <- grid.arrange(plot1, plot2, plot3, ncol = 1, heights = c(1,1,1))

library(cowplot)

# 提取plot2的图例
legend <- get_legend(plot0)

# 将图例添加到右侧
final_plot <- plot_grid(plot, legend, ncol = 2, rel_widths = c(0.9, 0.1))

# 显示最终的组合图
final_plot


ggsave("~/Desktop/广东省肿瘤数据/Plot(5)/按exposure区分的各癌种case条图/inc_Both_kind.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 12.7, 
       height = 7.2, 
       units = "in", 
       dpi = 1000)

#============================按exposure区分的各癌种case条图mort（更新）==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)
library(gridExtra)
data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019(1)/attribute_all_age.xlsx")

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 0)

# 定义癌种的顺序
site_order <- c(
  "Traches,bronchus and lung","Nasopharynx","Oral cavity, pharynx","Larynx",
  "Liver","Colorectum","Stomach","Esophagus","Pancreas","Gallbladder",
  "Cervix uteri","Breast","Corpus uteri","Prostate","Ovary","Vulva, Vagina", "Penis",
  "Thyroid","Bladder","Kidney","Myeloid leukaemia","Hodgkin disease"
)

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
data <- data %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序

# 定义初始的颜色节点
nodes <- c( "#bc6356","#dc917b","#eabaa1","#fee3ce")

# 创建插值函数
color_interpolation <- colorRampPalette(nodes)

# 生成22种颜色
colors1 <- color_interpolation(7)
print(colors1)

#custom_colors <- c("#3581b7","#4A92C3","#5FA3CF","#75B5DC",
"#9B3F5C" ,"#AD5873", "#BF718B", "#D08BA7", "#DFA9C5", "#EBCCE2",
"#BC6356" ,"#CC7968", "#DC917B", "#E3A58E", "#EABAA1", "#F3CEB7" ,"#FEE3CE",
"#4d7e54","#669877","#81b095","#a4cbb7","#cfeadf"
)
custom_colors <- c("#4FA8C7", "#6BB8D5", "#83C4E0", "#9ED8E7", 
                   "#D17D8E", "#D99CA8", "#E0A3B4", "#E7B8C0", "#F1C5D4", "#F5D0E1", 
                   "#D64F3A", "#E7725C", "#EB8A70", "#ED9D85", "#F1B1A3", "#F5C6B9","#F7D0C5",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8"
)
custom_colors <- c("#567CD0","#82A3E0","#A9C2EA" ,"#D6E8FA", 
                   "#F2A7C0", "#F5B6CC", "#F8C5D8", "#FBD4E4", "#FEE3ED", "#FFF2F7", 
                   "#E47A6D", "#E68E80", "#E8A292", "#EAB6A4", "#EDCAB6", "#EFD0A8", "#F1D6B9",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8"
)

data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

data <- data %>%
  mutate(sex = as.character(sex))

data <- data %>%
  group_by(exposure, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(mort_attribute))  #%>%


plot0 <- ggplot(data, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 18000), breaks = c(0, 6000, 12000, 18000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),  # 设置图例标题
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 设置图例文本样式
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -10000, xmax = -8000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示


data1 <- data %>% filter(sex == "Total")
exposure_order1 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Low vegetable intake","Helicobacter pylori","Second-hand smoking",
  "Human papillomavirus","Excess bodyweight","Red meat","Epstein-Barr virus",
  "Diabetes", "Low fruit intake", "Alcohol","Hepatitis B virus","Smoking"
)
data1 <- data1 %>%
  mutate(exposure = factor(exposure, levels = exposure_order1))  # 按照提供的顺序调整 site 的因子顺序

plot1 <- ggplot(data1, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 18000), breaks = c(0, 6000, 12000, 18000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = NULL, y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 去除x轴标签
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Total", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -4000, xmax = -3000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

data2 <- data %>% filter(sex == "Male")
exposure_order2 <- c(
  "Physical inactivity","Clonorchis sinensis","Hepatitis C virus","Low vegetable intake",
  "Human papillomavirus","Helicobacter pylori","Second-hand smoking","Excess bodyweight","Red meat","Epstein-Barr virus",
  "Diabetes", "Low fruit intake","Hepatitis B virus", "Alcohol","Smoking"
)
data2 <- data2 %>%
  mutate(exposure = factor(exposure, levels = exposure_order2))  # 按照提供的顺序调整 site 的因子顺序

plot2 <- ggplot(data2, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 18000), breaks = c(0, 6000, 12000, 18000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = "Exposure", fill = "Cancer Site") + 
  labs(y = NULL) +
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_blank(),  # 去除x轴标签
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Male", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -4000, xmax = -3000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data2, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

data3 <- data %>% filter(sex == "Female")
exposure_order3 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Low vegetable intake","Physical inactivity","Helicobacter pylori",
  "Smoking","Red meat","Second-hand smoking","Alcohol","Excess bodyweight","Epstein-Barr virus",
  "Hepatitis B virus","Diabetes","Human papillomavirus","Low fruit intake"
)
data3 <- data3 %>%
  mutate(exposure = factor(exposure, levels = exposure_order3))  # 按照提供的顺序调整 site 的因子顺序

plot3 <- ggplot(data3, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 18000), breaks = c(0, 6000, 12000, 18000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.text = element_blank(),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -4000, xmax = -3000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data3, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

plot1 <- plot1 + theme(plot.margin = margin(t = 10, r = 100, b = 1, l = 50))  # Adjust for plot1
plot2 <- plot2 + theme(plot.margin = margin(t = 5, r = 100, b = 1, l = 50))  # Adjust for plot2
plot3 <- plot3 + theme(plot.margin = margin(t = 5, r = 100, b = 10, l = 50))  # Increase bottom margin for plot3

# 使用 grid.arrange() 将三个图排列
plot <- grid.arrange(plot1, plot2, plot3, ncol = 1, heights = c(1,1,1))

library(cowplot)

# 提取plot2的图例
legend <- get_legend(plot0)

# 将图例添加到右侧
final_plot <- plot_grid(plot, legend, ncol = 2, rel_widths = c(0.9, 0.1))

# 显示最终的组合图
final_plot


ggsave("~/Desktop/广东省肿瘤数据/Plot(5)/按exposure区分的各癌种case条图/mort_Both_kind.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 12.7, 
       height = 7.2, 
       units = "in", 
       dpi = 1000)


#============================按exposure区分的各癌种case条图inc-城市（更新）==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)

data <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_all_age.xlsx")

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 1)

# 定义癌种的顺序
site_order <- c(
  "Traches,bronchus and lung","Nasopharynx","Oral cavity, pharynx","Larynx",
  "Liver","Colorectum","Stomach","Esophagus","Pancreas","Gallbladder",
  "Cervix uteri","Breast","Corpus uteri","Prostate","Ovary","Vulva, Vagina", "Penis",
  "Thyroid","Bladder","Kidney","Myeloid leukaemia","Hodgkin disease"
)

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
data <- data %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序

# 定义初始的颜色节点
nodes <- c( "#bc6356","#dc917b","#eabaa1","#fee3ce")

# 创建插值函数
color_interpolation <- colorRampPalette(nodes)

# 生成22种颜色
colors1 <- color_interpolation(7)
print(colors1)

custom_colors <- c("#3581b7","#4A92C3","#5FA3CF","#75B5DC",
                   "#9B3F5C" ,"#AD5873", "#BF718B", "#D08BA7", "#DFA9C5", "#EBCCE2",
                   "#BC6356" ,"#CC7968", "#DC917B", "#E3A58E", "#EABAA1", "#F3CEB7" ,"#FEE3CE",
                   "#4d7e54","#669877","#81b095","#a4cbb7","#cfeadf"
)

custom_colors <- c("#567CD0","#82A3E0","#A9C2EA" ,"#D6E8FA", 
                   "#F2A7C0", "#F5B6CC", "#F8C5D8", "#FBD4E4", "#FEE3ED", "#FFF2F7", 
                   "#E47A6D", "#E68E80", "#E8A292", "#EAB6A4", "#EDCAB6", "#EFD0A8", "#F1D6B9",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8"
)

data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

data <- data %>%
  mutate(sex = as.character(sex))

data <- data %>%
  group_by(exposure, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(inc_attribute))  #%>%


plot0 <- ggplot(data, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 12000), breaks = c(0,3000, 6000, 9000, 12000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),  # 设置图例标题
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 设置图例文本样式
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -2600, xmax = -2000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示



data1 <- data %>% filter(sex == "Total")
exposure_order1 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Helicobacter pylori","Second-hand smoking","Low vegetable intake",
  "Red meat","Excess bodyweight","Human papillomavirus","Epstein-Barr virus",
  "Low fruit intake","Hepatitis B virus", "Diabetes", "Alcohol","Smoking"
)
data1 <- data1 %>%
  mutate(exposure = factor(exposure, levels = exposure_order1))  # 按照提供的顺序调整 site 的因子顺序

plot1 <- ggplot(data1, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 12000), breaks = c(0,3000, 6000, 9000, 12000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = NULL, y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 去除x轴标签
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Total", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -2600, xmax = -2000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

data2 <- data %>% filter(sex == "Male")
exposure_order2 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Second-hand smoking",
  "Human papillomavirus","Helicobacter pylori","Low vegetable intake","Red meat","Excess bodyweight","Epstein-Barr virus",
  "Low fruit intake","Diabetes","Hepatitis B virus", "Alcohol","Smoking"
)
data2 <- data2 %>%
  mutate(exposure = factor(exposure, levels = exposure_order2))  # 按照提供的顺序调整 site 的因子顺序

plot2 <- ggplot(data2, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 12000), breaks = c(0,3000, 6000, 9000, 12000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = "Exposure", fill = "Cancer Site") + 
  labs(y = NULL) +
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_blank(),  # 去除x轴标签
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Male", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -2600, xmax = -2000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data2, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

data3 <- data %>% filter(sex == "Female")
exposure_order3 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Low vegetable intake","Helicobacter pylori",
  "Smoking","Red meat","Second-hand smoking","Hepatitis B virus","Epstein-Barr virus","Alcohol",
  "Low fruit intake", "Excess bodyweight","Diabetes","Human papillomavirus"
)
data3 <- data3 %>%
  mutate(exposure = factor(exposure, levels = exposure_order3))  # 按照提供的顺序调整 site 的因子顺序

plot3 <- ggplot(data3, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 12000), breaks = c(0,3000, 6000, 9000, 12000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.text = element_blank(),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -2600, xmax = -2000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data3, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

plot1 <- plot1 + theme(plot.margin = margin(t = 10, r = 100, b = 1, l = 50))  # Adjust for plot1
plot2 <- plot2 + theme(plot.margin = margin(t = 5, r = 100, b = 1, l = 50))  # Adjust for plot2
plot3 <- plot3 + theme(plot.margin = margin(t = 5, r = 100, b = 10, l = 50))  # Increase bottom margin for plot3

# 使用 grid.arrange() 将三个图排列
plot <- grid.arrange(plot1, plot2, plot3, ncol = 1, heights = c(1,1,1))

library(cowplot)

# 提取plot2的图例
legend <- get_legend(plot0)

# 将图例添加到右侧
final_plot <- plot_grid(plot, legend, ncol = 2, rel_widths = c(0.9, 0.1))

# 显示最终的组合图
final_plot


ggsave("~/Desktop/广东省肿瘤数据/Plot(5)/按exposure区分的各癌种case条图/inc_Urban.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 12.7, 
       height = 7.2, 
       units = "in", 
       dpi = 1000)

#============================按exposure区分的各癌种case条图mort-城市（更新）==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)

data <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_all_age.xlsx")

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 1)

# 定义癌种的顺序
site_order <- c(
  "Traches,bronchus and lung","Nasopharynx","Oral cavity, pharynx","Larynx",
  "Liver","Colorectum","Stomach","Esophagus","Pancreas","Gallbladder",
  "Cervix uteri","Breast","Corpus uteri","Prostate","Ovary","Vulva, Vagina", "Penis",
  "Thyroid","Bladder","Kidney","Myeloid leukaemia","Hodgkin disease"
)

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
data <- data %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序

# 定义初始的颜色节点
nodes <- c( "#bc6356","#dc917b","#eabaa1","#fee3ce")

# 创建插值函数
color_interpolation <- colorRampPalette(nodes)

# 生成22种颜色
colors1 <- color_interpolation(7)
print(colors1)

custom_colors <- c("#3581b7","#4A92C3","#5FA3CF","#75B5DC",
                   "#9B3F5C" ,"#AD5873", "#BF718B", "#D08BA7", "#DFA9C5", "#EBCCE2",
                   "#BC6356" ,"#CC7968", "#DC917B", "#E3A58E", "#EABAA1", "#F3CEB7" ,"#FEE3CE",
                   "#4d7e54","#669877","#81b095","#a4cbb7","#cfeadf"
)

custom_colors <- c("#567CD0","#82A3E0","#A9C2EA" ,"#D6E8FA", 
                   "#F2A7C0", "#F5B6CC", "#F8C5D8", "#FBD4E4", "#FEE3ED", "#FFF2F7", 
                   "#E47A6D", "#E68E80", "#E8A292", "#EAB6A4", "#EDCAB6", "#EFD0A8", "#F1D6B9",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8"
)

data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

data <- data %>%
  mutate(sex = as.character(sex))

data <- data %>%
  group_by(exposure, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(mort_attribute))  #%>%


plot0 <- ggplot(data, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 9000), breaks = c(0,3000, 6000, 9000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),  # 设置图例标题
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 设置图例文本样式
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -2600, xmax = -2000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示



data1 <- data %>% filter(sex == "Total")
exposure_order1 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity",
  "Low vegetable intake","Helicobacter pylori","Second-hand smoking",
  "Human papillomavirus",
  "Red meat","Excess bodyweight","Epstein-Barr virus",
  "Low fruit intake", "Diabetes","Hepatitis B virus", "Alcohol","Smoking"
)
data1 <- data1 %>%
  mutate(exposure = factor(exposure, levels = exposure_order1))  # 按照提供的顺序调整 site 的因子顺序

plot1 <- ggplot(data1, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 9000), breaks = c(0,3000, 6000, 9000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = NULL, y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 去除x轴标签
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Total", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -1800, xmax = -1500, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

data2 <- data %>% filter(sex == "Male")
exposure_order2 <- c(
  "Physical inactivity","Clonorchis sinensis","Hepatitis C virus","Human papillomavirus","Second-hand smoking",
  "Helicobacter pylori","Low vegetable intake","Excess bodyweight","Red meat","Epstein-Barr virus",
  "Low fruit intake","Diabetes","Hepatitis B virus", "Alcohol","Smoking"
)
data2 <- data2 %>%
  mutate(exposure = factor(exposure, levels = exposure_order2))  # 按照提供的顺序调整 site 的因子顺序

plot2 <- ggplot(data2, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 9000), breaks = c(0,3000, 6000, 9000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = "Exposure", fill = "Cancer Site") + 
  labs(y = NULL) +
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_blank(),  # 去除x轴标签
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Male", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -1800, xmax = -1500, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data2, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

data3 <- data %>% filter(sex == "Female")
exposure_order3 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Low vegetable intake","Helicobacter pylori",
  "Smoking","Red meat","Second-hand smoking","Epstein-Barr virus","Alcohol",
  "Excess bodyweight","Hepatitis B virus","Low fruit intake","Human papillomavirus", "Diabetes"
)
data3 <- data3 %>%
  mutate(exposure = factor(exposure, levels = exposure_order3))  # 按照提供的顺序调整 site 的因子顺序

plot3 <- ggplot(data3, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 9000), breaks = c(0,3000, 6000, 9000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.text = element_blank(),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -1800, xmax = -1500, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data3, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

plot1 <- plot1 + theme(plot.margin = margin(t = 10, r = 100, b = 1, l = 50))  # Adjust for plot1
plot2 <- plot2 + theme(plot.margin = margin(t = 5, r = 100, b = 1, l = 50))  # Adjust for plot2
plot3 <- plot3 + theme(plot.margin = margin(t = 5, r = 100, b = 10, l = 50))  # Increase bottom margin for plot3

# 使用 grid.arrange() 将三个图排列
plot <- grid.arrange(plot1, plot2, plot3, ncol = 1, heights = c(1,1,1))

library(cowplot)

# 提取plot2的图例
legend <- get_legend(plot0)

# 将图例添加到右侧
final_plot <- plot_grid(plot, legend, ncol = 2, rel_widths = c(0.9, 0.1))

# 显示最终的组合图
final_plot


ggsave("~/Desktop/广东省肿瘤数据/Plot(5)/按exposure区分的各癌种case条图/mort_Urban.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 12.7, 
       height = 7.2, 
       units = "in", 
       dpi = 1000)
#============================按exposure区分的各癌种case条图inc-农村（更新）==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)

data <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_all_age.xlsx")

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 2)

# 定义癌种的顺序
site_order <- c(
  "Traches,bronchus and lung","Nasopharynx","Oral cavity, pharynx","Larynx",
  "Liver","Colorectum","Stomach","Esophagus","Pancreas","Gallbladder",
  "Cervix uteri","Breast","Corpus uteri","Prostate","Ovary","Vulva, Vagina", "Penis",
  "Thyroid","Bladder","Kidney","Myeloid leukaemia","Hodgkin disease"
)

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
data <- data %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序

# 定义初始的颜色节点
nodes <- c( "#bc6356","#dc917b","#eabaa1","#fee3ce")

# 创建插值函数
color_interpolation <- colorRampPalette(nodes)

# 生成22种颜色
colors1 <- color_interpolation(7)
print(colors1)

custom_colors <- c("#3581b7","#4A92C3","#5FA3CF","#75B5DC",
                   "#9B3F5C" ,"#AD5873", "#BF718B", "#D08BA7", "#DFA9C5", "#EBCCE2",
                   "#BC6356" ,"#CC7968", "#DC917B", "#E3A58E", "#EABAA1", "#F3CEB7" ,"#FEE3CE",
                   "#4d7e54","#669877","#81b095","#a4cbb7","#cfeadf"
)

custom_colors <- c("#567CD0","#82A3E0","#A9C2EA" ,"#D6E8FA", 
                   "#F2A7C0", "#F5B6CC", "#F8C5D8", "#FBD4E4", "#FEE3ED", "#FFF2F7", 
                   "#E47A6D", "#E68E80", "#E8A292", "#EAB6A4", "#EDCAB6", "#EFD0A8", "#F1D6B9",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8"
)

data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

data <- data %>%
  mutate(sex = as.character(sex))

data <- data %>%
  group_by(exposure, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(inc_attribute))  #%>%


plot0 <- ggplot(data, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 12500), breaks = c(0,3000, 6000, 9000, 12000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),  # 设置图例标题
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 设置图例文本样式
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -2650, xmax = -2000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示



data1 <- data %>% filter(sex == "Total")
exposure_order1 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Low vegetable intake","Physical inactivity","Helicobacter pylori","Second-hand smoking",
  "Excess bodyweight","Red meat","Human papillomavirus","Epstein-Barr virus","Diabetes", 
  "Hepatitis B virus", "Low fruit intake","Alcohol","Smoking"
)
data1 <- data1 %>%
  mutate(exposure = factor(exposure, levels = exposure_order1))  # 按照提供的顺序调整 site 的因子顺序

plot1 <- ggplot(data1, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 12500), breaks = c(0,3000, 6000, 9000, 12000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = NULL, y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 去除x轴标签
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Total", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -2720, xmax = -2000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

data2 <- data %>% filter(sex == "Male")
exposure_order2 <- c(
  "Clonorchis sinensis","Low vegetable intake","Hepatitis C virus","Physical inactivity",
  "Human papillomavirus","Second-hand smoking","Helicobacter pylori","Excess bodyweight","Red meat","Diabetes","Epstein-Barr virus",
  "Low fruit intake","Hepatitis B virus", "Alcohol","Smoking"
)
data2 <- data2 %>%
  mutate(exposure = factor(exposure, levels = exposure_order2))  # 按照提供的顺序调整 site 的因子顺序

plot2 <- ggplot(data2, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 12500), breaks = c(0,3000, 6000, 9000, 12000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = "Exposure", fill = "Cancer Site") + 
  labs(y = NULL) +
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_blank(),  # 去除x轴标签
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Male", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -2720, xmax = -2000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data2, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

data3 <- data %>% filter(sex == "Female")
exposure_order3 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Low vegetable intake","Physical inactivity","Helicobacter pylori",
  "Smoking","Alcohol","Red meat","Second-hand smoking","Hepatitis B virus","Epstein-Barr virus",
  "Excess bodyweight","Diabetes","Low fruit intake","Human papillomavirus"
)
data3 <- data3 %>%
  mutate(exposure = factor(exposure, levels = exposure_order3))  # 按照提供的顺序调整 site 的因子顺序

plot3 <- ggplot(data3, aes(y = exposure, x = inc_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 12500), breaks = c(0,3000, 6000, 9000, 12000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.text = element_blank(),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -2720, xmax = -2000, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data3, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

plot1 <- plot1 + theme(plot.margin = margin(t = 10, r = 100, b = 1, l = 50))  # Adjust for plot1
plot2 <- plot2 + theme(plot.margin = margin(t = 5, r = 100, b = 1, l = 50))  # Adjust for plot2
plot3 <- plot3 + theme(plot.margin = margin(t = 5, r = 100, b = 10, l = 50))  # Increase bottom margin for plot3

# 使用 grid.arrange() 将三个图排列
plot <- grid.arrange(plot1, plot2, plot3, ncol = 1, heights = c(1,1,1))

library(cowplot)

# 提取plot2的图例
legend <- get_legend(plot0)

# 将图例添加到右侧
final_plot <- plot_grid(plot, legend, ncol = 2, rel_widths = c(0.9, 0.1))

# 显示最终的组合图
final_plot

ggsave("~/Desktop/广东省肿瘤数据/Plot(5)/按exposure区分的各癌种case条图/inc_Rural.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 12.7, 
       height = 7.2, 
       units = "in", 
       dpi = 1000)

#============================按exposure区分的各癌种case条图mort-农村（更新）==============================================
install.packages("readxl") 
install.packages("ggplot2")
install.packages("dplyr")

library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)

data <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_all_age.xlsx")

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 2)

# 定义癌种的顺序
site_order <- c(
  "Traches,bronchus and lung","Nasopharynx","Oral cavity, pharynx","Larynx",
  "Liver","Colorectum","Stomach","Esophagus","Pancreas","Gallbladder",
  "Cervix uteri","Breast","Corpus uteri","Prostate","Ovary","Vulva, Vagina", "Penis",
  "Thyroid","Bladder","Kidney","Myeloid leukaemia","Hodgkin disease"
)

# 在绘制堆积条形图之前，确保 site 按照指定的顺序排序
data <- data %>%
  mutate(site = factor(site, levels = site_order))  # 按照提供的顺序调整 site 的因子顺序

# 定义初始的颜色节点
nodes <- c( "#bc6356","#dc917b","#eabaa1","#fee3ce")

# 创建插值函数
color_interpolation <- colorRampPalette(nodes)

# 生成22种颜色
colors1 <- color_interpolation(7)
print(colors1)

custom_colors <- c("#3581b7","#4A92C3","#5FA3CF","#75B5DC",
                   "#9B3F5C" ,"#AD5873", "#BF718B", "#D08BA7", "#DFA9C5", "#EBCCE2",
                   "#BC6356" ,"#CC7968", "#DC917B", "#E3A58E", "#EABAA1", "#F3CEB7" ,"#FEE3CE",
                   "#4d7e54","#669877","#81b095","#a4cbb7","#cfeadf"
)

custom_colors <- c("#567CD0","#82A3E0","#A9C2EA" ,"#D6E8FA", 
                   "#F2A7C0", "#F5B6CC", "#F8C5D8", "#FBD4E4", "#FEE3ED", "#FFF2F7", 
                   "#E47A6D", "#E68E80", "#E8A292", "#EAB6A4", "#EDCAB6", "#EFD0A8", "#F1D6B9",
                   "#6ABF7B", "#82D6A0", "#A2E0A5", "#C6E7B9", "#E3F5E8"
)


data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

data <- data %>%
  mutate(sex = as.character(sex))

data <- data %>%
  group_by(exposure, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(mort_attribute))  #%>%


plot0 <- ggplot(data, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 12500), breaks = c(0,3000, 6000, 9000, 12000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),  # 设置图例标题
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 设置图例文本样式
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -1800, xmax = -1500, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示



data1 <- data %>% filter(sex == "Total")
exposure_order1 <- c(
  "Low vegetable intake","Clonorchis sinensis","Hepatitis C virus","Physical inactivity","Helicobacter pylori","Excess bodyweight","Second-hand smoking",
  "Human papillomavirus","Red meat","Epstein-Barr virus","Diabetes", 
  "Low fruit intake","Alcohol","Hepatitis B virus","Smoking"
)
data1 <- data1 %>%
  mutate(exposure = factor(exposure, levels = exposure_order1))  # 按照提供的顺序调整 site 的因子顺序

plot1 <- ggplot(data1, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 9000), breaks = c(0,3000, 6000, 9000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = NULL, y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 去除x轴标签
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Total", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -1800, xmax = -1500, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data1, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

data2 <- data %>% filter(sex == "Male")
exposure_order2 <- c(
  "Low vegetable intake","Physical inactivity","Clonorchis sinensis","Hepatitis C virus",
  "Human papillomavirus","Excess bodyweight","Helicobacter pylori","Second-hand smoking","Red meat","Epstein-Barr virus",
  "Diabetes","Low fruit intake","Alcohol","Hepatitis B virus", "Smoking"
)
data2 <- data2 %>%
  mutate(exposure = factor(exposure, levels = exposure_order2))  # 按照提供的顺序调整 site 的因子顺序

plot2 <- ggplot(data2, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 9000), breaks = c(0,3000, 6000, 9000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = "Exposure", fill = "Cancer Site") + 
  labs(y = NULL) +
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_blank(),  # 去除x轴标签
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_blank(),  # 去除x轴标题
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 8, family = "Times New Roman", color = "black"),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签 
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Male", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -1800, xmax = -1500, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data2, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

data3 <- data %>% filter(sex == "Female")
exposure_order3 <- c(
  "Clonorchis sinensis","Hepatitis C virus","Low vegetable intake","Physical inactivity","Helicobacter pylori",
  "Excess bodyweight","Alcohol","Smoking","Red meat","Second-hand smoking","Epstein-Barr virus",
  "Diabetes","Human papillomavirus","Hepatitis B virus","Low fruit intake"
)
data3 <- data3 %>%
  mutate(exposure = factor(exposure, levels = exposure_order3))  # 按照提供的顺序调整 site 的因子顺序

plot3 <- ggplot(data3, aes(y = exposure, x = mort_attribute, fill = site)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 9000), breaks = c(0,3000, 6000, 9000), expand = c(0, 0)) +  # 扩展x轴范围
  labs(x = "Number of Attributable Cases", y = NULL, fill = "Cancer Site") + 
  scale_fill_manual(values = custom_colors, guide = guide_legend(ncol = 1)) +  # 使用自定义颜色
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 8, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.title = element_blank(),  # 删除图例标题
        legend.text = element_blank(),  # 图例文本样式
        legend.position = "none",  # 隐藏图例
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines"),
        plot.margin = margin(l = 50)) + # 调整左边距以腾出空间给文本标签
  # 使用annotation_custom来添加纵坐标标签左侧的文本
  annotation_custom(
    grob = textGrob("Female", gp = gpar(fontsize = 11, fontfamily = "Times New Roman", col = "black", fontface = "bold"),rot = 0),  # 添加文本并使用 gpar 设置样式
    xmin = -1800, xmax = -1500, ymin = 0, ymax = 16  # 设置文本的坐标
  )+
  # 在条形图右侧添加总和标签
  geom_text(data = data3, aes(x = total_cases+50, y = exposure, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示

plot1 <- plot1 + theme(plot.margin = margin(t = 10, r = 100, b = 1, l = 50))  # Adjust for plot1
plot2 <- plot2 + theme(plot.margin = margin(t = 5, r = 100, b = 1, l = 50))  # Adjust for plot2
plot3 <- plot3 + theme(plot.margin = margin(t = 5, r = 100, b = 10, l = 50))  # Increase bottom margin for plot3

# 使用 grid.arrange() 将三个图排列
plot <- grid.arrange(plot1, plot2, plot3, ncol = 1, heights = c(1,1,1))

library(cowplot)

# 提取plot2的图例
legend <- get_legend(plot0)

# 将图例添加到右侧
final_plot <- plot_grid(plot, legend, ncol = 2, rel_widths = c(0.9, 0.1))

# 显示最终的组合图
final_plot

ggsave("~/Desktop/广东省肿瘤数据/Plot(5)/按exposure区分的各癌种case条图/mort_Rural.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 12.7, 
       height = 7.2, 
       units = "in", 
       dpi = 1000)
#===================按factor区分的各癌种PAF条图(合并Total、Urban、Rural)==============================================

data <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/各个factor分癌种+所有癌症的PAF.xlsx")
data1 <- read_excel("~/Desktop/广东省肿瘤数据/PAF2019/各个factor分癌种的PAF.xlsx")
data2 <- read_xlsx("~/Desktop/广东省肿瘤数据/PAF2019/attribute_case_all_exposure_ICD10(all_age).xlsx")
data2 <- data2 %>%
  mutate(factor = "All")

data <- rbind(data1, data2)


panel_labels_factor <- c("Infectious agent" = "Infectious agent",
                         "Metabolic factor" = "Metabolic factor",
                         "Dietary factor" = "Dietary factor", 
                         "Behaviour factor" = "Behaviour factor", 
                         "All" = "All factor")

panel_labels_kind <- c("2" = "Rural", 
                       "1" = "Urban", 
                       "0" = "Total")


panel_labels_sex <- c("2" = "Female", 
                      "1" = "Male", 
                      "0" = "Total")

data$PAF_inc <- as.numeric(data$PAF_inc)
data$PAF_inc<-data$PAF_inc*100

data$mort_attribute <- as.integer(data$mort_attribute)
data$mort_attribute_lower <- as.integer(data$mort_attribute_lower)
data$mort_attribute_upper <- as.integer(data$mort_attribute_upper)
data$inc_attribute <- as.integer(data$inc_attribute)
data$inc_attribute_lower <- as.integer(data$inc_attribute_lower)
data$inc_attribute_upper <- as.integer(data$inc_attribute_upper)

data <- data %>% filter(kind == 0)


data <- data %>%
  group_by(factor, sex) %>%  # 按exposure和sex分组
  mutate(total_cases = sum(inc_attribute))  #%>%

site_colors <- viridis(n = length(unique(data$site)))  # 根据site列的唯一值数量生成颜色

ggplot(data, aes(y = factor, x = PAF_inc, fill = site)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +  # 堆积条形图
  scale_x_continuous(name = "PAF", 
                     limits = c(0, 100),  # 设置x轴范围
                     labels = scales::comma_format(),  # 去掉百分号，改为普通数字格式
                     expand = c(0, 0)) +  # 横坐标范围
  scale_fill_manual(values = site_colors) +  # 使用viridis调色板的颜色为site分配颜色
  labs(x = "PAF", y = NULL) + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 10, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 10, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.position = "right",  # 图例放置在右侧
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 10, family = "Times New Roman", color = "black"),  # 图例文本样式
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside",  # 确保面板标签在外部
        strip.background = element_blank(),  # 清除面板标签的背景色
        panel.spacing.x = unit(0.5, "lines"), # 调整面板间距
        panel.spacing.y = unit(0.5, "lines")) +
  facet_wrap(~ sex, ncol = 1, scales = "free_y", 
             labeller = as_labeller(panel_labels_sex), 
             strip.position = "left")+  # 将面板标签放置在左侧
  guides(fill = guide_legend(ncol = 1))+  # 将图例设置为一列
  geom_text(data = data, aes(x = 60, y = factor, label = total_cases), 
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")+
  coord_cartesian(clip = "off")  # 允许超出范围的部分显示



#===================最新-按factor区分的PAF条图(合并Total、Urban、Rural)==============================================
library(readxl) 
library(ggplot2)
library(dplyr)
library(viridis)
library(scales)

data <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/某一factor所有癌症PAF及case.xlsx")
data1 <- read_xlsx("~/Desktop/数据/广东省肿瘤归因/分析结果/PAF2019/PAF2019(重新合并男女城乡)/所有危险因素-所有癌症PAF及case.xlsx")
data1 <- data1 %>%
  mutate(factor = "All")
data<-rbind(data,data1)

#data <- data %>% filter(kind == 0)

data$PAF_inc<-data$PAF_inc*100
data$PAF_mort<-data$PAF_mort*100

# 定义暴露因素的顺序
factor_order <- c("All","Infectious agent", "Metabolic factor","Dietary factor","Behaviour factor")

# 确保每个癌种都包含所有暴露因素，并按暴露因素的顺序排列
data <- data %>%
  mutate(factor = factor(factor, levels = factor_order)
  )  # 确保暴露因子按顺序排列

data$sex <- factor(data$sex, levels = c(2, 1, 0), labels = c("Female", "Male","Total" ))

#"Male" = "#7fb2d5", "Female" = "#f47f72" , "Total" = "#f6cd96""#fbd178""#abd9d0""#ABADAF""#f6cd96"

panel_labels <- c("0" = "Total", "1" = "Urban", "2" = "Rural")

ggplot(data, aes(y = factor, x = PAF_mort, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +  # 绘制条形图
  scale_x_continuous(name = "PAF(%)", 
                     limits = c(0, 60),
                     breaks = c(0, 20, 40, 60), 
                     labels = scales::comma_format(),  # 去掉百分号，改为普通数字格式
                     expand = c(0, 0)) +  # 横坐标范围
  scale_fill_manual(values = c("Male" = "#3366CC", "Female" = "#FF99CC" , "Total" = "#D0D4D7")) +
  labs(x = "PAF(%)", y = NULL) + 
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 10, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 10, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.position = "right",  # 图例放置在右侧
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.3, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 9, family = "Times New Roman", color = "black"),  # 图例文本样式
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank(),
        strip.text = element_text(size = 9, family = "Times New Roman", color = "black", face = "bold"),  # 设置facet标签的样式
        strip.placement = "outside") +  # 设置面板标签在纵坐标轴外侧
  geom_text(aes(x = PAF_mort + 1, y = factor, label = sprintf("%.2f", PAF_mort)),  # 保留一位小数
            position = position_dodge(width = 0.8),  # 保证标签的位置与条形图对齐
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black") +
  # 在x = 100的位置添加inc_95CI文本标签
  geom_text(aes(x = 105, y = factor, label = mort_95CI),  # 设置x = 100，确保文本标签在横坐标轴外侧
            position = position_dodge(width = 0.8),  # 确保文本与条形图对齐
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black") +
  facet_wrap(~ kind, ncol = 1, scales = "free_y", 
             labeller = as_labeller(panel_labels), 
             strip.position = "left")  # 将面板标签放置在左侧

ggsave("~/Desktop/广东省肿瘤数据/Plot(5)/按factor区分的各癌种case条图/mort-Both_kind.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 8, 
       height = 7, 
       units = "in", 
       dpi = 1000)

#=======================按factor区分的PAF条图==============================================
ggplot(data, aes(y = factor, x = PAF_inc, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +  # 绘制条形图
  scale_x_continuous(name = "PAF", 
                     limits = c(0, 100), 
                     labels = scales::percent_format(scale = 1), 
                     expand = c(0, 0)) +  # 横坐标范围
  scale_fill_manual(values = c("Male" = "#3F62A4", "Female" = "#AB1F2A", "Total" = "#707070")) +
  labs(x = "PAF", y = "factor") + 
  labs(y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
        axis.text.x = element_text(size = 10, family = "Times New Roman", color = "black"),
        axis.line.x = element_line(color = "grey", size = 0.3),  # 添加横坐标轴线
        axis.line.y = element_line(color = "grey", size = 0.3),  # 添加纵坐标线
        axis.title.x = element_text(size = 10, family = "Times New Roman", face = "bold", color = "black"),
        axis.ticks = element_line(color = "grey", size = 0.3),  # 添加坐标轴刻度线
        legend.position = "right",  # 图例放置在右侧
        legend.title = element_blank(),  # 删除图例标题
        legend.key.size = unit(0.4, "cm"),  # 缩小图例的大小
        legend.text = element_text(size = 10, family = "Times New Roman", color = "black"),  # 图例文本样式
        panel.grid = element_blank(),  
        panel.background = element_blank(),  
        plot.background = element_blank()) +
  geom_text(data = data, 
            aes(x = PAF_inc + 1, y = factor, label = sprintf("%.1f", PAF_inc)),  # 保留两位小数
            position = position_dodge(width = 0.8),  # 保证标签的位置与条形图对齐
            vjust = 0.5, hjust = 0, size = 3, family = "Times New Roman", color = "black")

ggsave("~/Desktop/广东省肿瘤数据/Plot(4)/按factor区分的各癌种case条图/Both_kind.tiff", 
       plot = last_plot(), 
       device = "tiff", 
       width = 7, 
       height = 7, 
       units = "in", 
       dpi = 1000)

#============================================================================================================
data$PAF_inc<-data$PAF_inc*100
data$PAF_lower_inc<-data$PAF_lower_inc*100
data$PAF_upper_inc<-data$PAF_upper_inc*100

data$sex <- factor(data$sex, levels = c(0, 1, 2), labels = c("Total", "Male", "Female"))

data <- data %>%
  filter(
    (site %in% c("Vulva, Vagina", "Breast", "Ovary", "Corpus uteri", "Cervix uteri") & sex == "Female") |
      (site %in% c("Penis", "Prostate") & sex == "Male") |
      !(site %in% c("Vulva, Vagina", "Breast", "Ovary", "Corpus uteri", "Cervix uteri", "Penis", "Prostate"))
  )

# 定义癌种的顺序
site_order <- c(
  "Penis", "Vulva, Vagina", "Prostate", "Breast", "Ovary", "Corpus uteri", "Cervix uteri",  "Non-hodgkin lymphoma", "Hodgkin disease","Myeloid leukaemia", "Thyroid","Bladder", "Kidney", "Gallbladder", "Liver", "Pancreas", "Rectum", "Colon", "Stomach","Esophagus", "Traches,bronchus and lung", "Larynx", "Nasopharynx", "Oral cavity, pharynx"
)

ggplot() +
  # 绘制窄条形图（包括所有性别）
  geom_bar(data = narrow_data, aes(y = site, x = PAF_inc, fill = sex), 
           stat = "identity", position = position_dodge(width = 0), width = 0.25, show.legend = TRUE) +
  
  # 绘制宽条形图（包括所有性别）
  geom_bar(data = wide_data, aes(y = site, x = PAF_inc, fill = sex), 
           stat = "identity", position = position_dodge(width = 0.8), width = 0.8, show.legend = TRUE) +
  
  # 添加误差条（宽条形图部分）
  geom_errorbar(data = wide_data, aes(x = PAF_inc, y = site, xmin = PAF_lower_inc, xmax = PAF_upper_inc, group = sex), 
                position = position_dodge(width = 0.8), width = 0.3, color = "black") +
  
  # 添加误差条（窄条形图部分）
  geom_errorbar(data = narrow_data, aes(x = PAF_inc, y = site, xmin = PAF_lower_inc, xmax = PAF_upper_inc, group = sex), 
                position = position_dodge(width = 0.1), width = 0.1, color = "black") +
  
  # 设置 x 轴
  scale_x_continuous(name = "PAF", limits = c(0, 100), 
                     labels = scales::percent_format(scale = 1), expand = c(0, 0)) +
  
  # 设置填充颜色
  scale_fill_manual(values = c("Female" = "#CD0000", "Male" = "#1874CD", "Total" = "#8B8989")) +
  
  # 不设置图例标题
  labs(y = NULL, fill = NULL) +
  
  # 使用简洁主题
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 9, family = "Times New Roman", color = "black"),
    axis.text.y = element_text(size = 9, family = "Times New Roman", color = "black", angle = 0, hjust = 1),  
    axis.line.x = element_line(color = "grey", size = 0.3), 
    axis.line.y = element_line(color = "grey", size = 0.3), 
    axis.title.x = element_text(size = 9, family = "Times New Roman", face = "bold", color = "black"),  
    axis.ticks = element_line(color = "grey", size = 0.3),
    legend.position = "right",  # 设置图例位置
    legend.title = element_blank(),  # 删除图例标题
    legend.key.size = unit(0.4, "cm"),  # 缩小图例的大小
    legend.text = element_text(size = 9, family = "Times New Roman", color = "black"),  # 图例文本样式
    panel.grid = element_blank(),  
    panel.background = element_blank(),  
    plot.background = element_blank()
  )


#=========================================
#NPC_EB
meta <- read_excel("~/Desktop/广东省肿瘤数据/暴露/感染因素-暴露.xlsx", sheet = "EB-NPC")
meta <- meta[,c(1,7,8,10,11)]

NPC_EB <- read_excel("~/Desktop/HPV/excel/NPC&EB.xlsx", sheet = "广东1990-2017")
NPC_EB <-  NPC_EB[,c(1,2)]
NPC_EB <- NPC_EB %>%
  distinct(Study, .keep_all = TRUE)

meta<-meta%>%
  left_join(NPC_EB,by="Study")

NPC_EB <- read_excel("~/Desktop/HPV/excel/NPC&EB.xlsx", sheet = "广东2018-2024")
NPC_EB <-  NPC_EB[,c(1,2)]
NPC_EB <- NPC_EB %>%
  distinct(Study, .keep_all = TRUE)

meta<-meta%>%
  left_join(NPC_EB,by="Study")

library(openxlsx)
write.xlsx(meta,"~/Desktop/广东省肿瘤数据/暴露/meta_NPC_EB.xlsx")





library(dplyr)
library(tidyr)
meta <- meta %>%
  unite(ID, id, id2, sep = "", remove = TRUE, na.rm = TRUE)

install.packages("colorspace")

# 加载 colorspace 包
library(colorspace)

"#478ecc""#F26B8F"

color <-"#F26B8F"
lighter_color <- lighten(color, amount = 0.3)  # 亮度增加
lighter_color



