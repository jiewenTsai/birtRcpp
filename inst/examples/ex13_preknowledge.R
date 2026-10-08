# Example 13. Detecting item preknowledge with response times (Sinharay, 2020, Table 2), on the licensure
# data of LNIRT (CredentialForm1, operational items 1-170). Test takers with a zero response time are
# left out, which gives exactly the 1,624 test takers and 41 flagged test takers of Sinharay (2020).
# Compromised items: the 64 items flagged in the public version of the data
# (github.com/czopluoglu/emip-mglnrt, dataset2); Sinharay used 63.
# Item parameters: ecm() of rtirt_null() (2PL scores and lognormal times jointly; Sinharay estimated the
# 2PL with mirt and the lognormal model with lavaan separately).
library(birtRcpp)
if (!requireNamespace("LNIRT", quietly = TRUE)) stop("this example uses the data of the LNIRT package")
data("CredentialForm1", package = "LNIRT")
cf <- get("CredentialForm1")
RT <- as.matrix(cf[paste0("idur.", 1:170)]); keep <- rowSums(RT == 0) == 0
cf <- cf[keep, ]
dat <- input_data(resp = as.matrix(cf[paste0("iraw.", 1:170)]), time = as.matrix(cf[paste0("idur.", 1:170)]))
compromised <- c(1, 3, 4, 7, 8, 9, 10, 14, 15, 21, 22, 25, 26, 29, 30, 35, 38, 40, 41, 43, 50, 53, 54, 58, 59, 61, 62, 67, 68, 72,
                 74, 75, 81, 84, 87, 89, 96, 97, 98, 106, 110, 116, 121, 127, 128, 129, 134, 135, 136, 140, 144, 145, 151, 152,
                 154, 159, 160, 161, 162, 163, 164, 165, 169, 170)
cat(sprintf("%d test takers, %d flagged; %d compromised items\n", nrow(cf), sum(cf$Flagged), length(compromised)))
t0 <- Sys.time(); fit <- ecm(rtirt_null(dat), se = FALSE); cat(sprintf("ecm(): %.0f s\n", as.numeric(Sys.time() - t0, units = "secs")))
pk <- preknowledge_test(fit, compromised)
pct <- function(p, g, alpha) round(100 * mean(p[g] < alpha))
fl <- cf$Flagged == 1
tab <- t(sapply(c(chi_pf = "p_chi_pf", Lambda_s = "p_Lambda_s", L_s = "p_L_s"), function(v)
  c(not_0.1 = pct(pk[[v]], !fl, .001), not_1 = pct(pk[[v]], !fl, .01), flagged_0.1 = pct(pk[[v]], fl, .001), flagged_1 = pct(pk[[v]], fl, .01))))
cat("Percent significant, Form 1 (not flagged at 0.1% / 1%, flagged at 0.1% / 1%)\n")
print(cbind(tab, sinharay = c("10 15 27 29", "1 3 22 29", "0 3 10 19")), quote = FALSE)
cat(sprintf("correlation of L_s and Lambda_s among the flagged: %.2f (Sinharay: .48)\n", cor(pk$L_s[fl], pk$Lambda_s[fl])))
