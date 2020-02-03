# # Create a free SalesForce account: https://developer.salesforce.com/signup

library(devtools)
library(salesforcer)

#we dont need to create the dict_to_df function because our outputs are already in df form. 

session <- function() {
  sf_auth(
    username = "your_username", 
    password = "your_password",
    security_token = "your_security_token"
  ) }

session <- session()

getAllFields <- function(objectName) {
  description <- rforcecom.getObjectDescription(session, objectName)
  fields <- as.character(description$name)
  rforcecom.retrieve(session, objectName, fields)
}

get_opportunities <- function(){
  objectName <- "Opportunity"
  fields <- c('CreatedDate','CloseDate','Name', 'StageName', 
              'ExpectedRevenue', 'Amount', 'LeadSource', 'IsWon', 
              'IsClosed', 'Type', 'Probability')
  rforcecom.retrieve(session, objectName, fields)}

get_cases <- function(){
  objectName <- "Case"
  fields <- c('CreatedDate', 'ContactId','Description', 
              'Type', 'Reason', 'Status', 'Origin', 'Subject', 
              'Priority', 'IsClosed', 'OwnerId', 'IsDeleted', 'AccountId')
  return(rforcecom.retrieve(session, objectName, fields))}

get_contacts <- function(){
  sf_Manager <- setClass('sf_Manager')
  objectName <- "Contact"
  fields <- c('Id', 'Salutation', 'FirstName', 'LastName')
  rforcecom.retrieve(session, objectName, fields)}

get_users <- function(){
  objectName <- "User"
  fields <- c('Id', 'FirstName', 'LastName')
  return(rforcecom.retrieve(session, objectName, fields))}

get_accounts <- function(){
  objectName <- "Account"
  fields <- c('Id', 'Name')
  return(rforcecom.retrieve(session, objectName, fields))}

add_opportunities <- function(name, stage, amount, probability, date, 
                              o_type, source){
  objectName <- "Opportunity"
  fields <- c(Name = name, StageName = stage, Amount = amount, 
              Probability = probability, CloseDate = date,
              Type = o_type, LeadSource = source)
  rforcecom.create(session, objectName, fields)}

add_leads <- function(company,status,state,source){
  objectName <- "Lead"
  fields <- c(LastName = company, Company = company, Status = status, 
              State = state, LeadSource = source)
  rforcecom.create(session, objectName, fields)}

add_cases <- function(account_id,origin,reason,subject,contact_id,case_type,
                      status,description,priority){
  objectName <- "Case"
  fields <- c(AccountId = account_id, Origin = origin, 
              Reason = reason, Subject = subject,
              ContactId = contact_id, Type = case_type, Status = status, 
              Description = description, Priority = priority)
  rforcecom.create(session, objectName, fields)}
