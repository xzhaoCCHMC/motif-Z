import smtplib
from email.mime.text import MIMEText as text

def send_email(sender,receiver,message,subject):
    sender = sender
    receivers = receiver
    m = text(message)
    m['Subject'] = subject
    m['From'] = sender
    m['To'] = receiver

   # message = message

    try:
       smtpObj = smtplib.SMTP('outbound-mail.cchmc.org')
       smtpObj.sendmail(sender, receivers, str(m))
       print "Successfully sent email"
    except SMTPException:
       print "Error: unable to send email" 

send_email('xueheng.zhao@cchmc.org', 'xueheng.zhao@gmail.com', 'hello', 'world')

