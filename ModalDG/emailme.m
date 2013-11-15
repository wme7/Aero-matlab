function [] = emailme(x)
%% Send email message when simulation finish to my Gmail account
% Based on:
% http://obasic.net/how-to-send-e-mail-from-matlab
%
% Coded by: Manuel Diaz

%% Set Parameters of my Gmail Account
setpref('Internet', 'E_mail', 'manuel.ade@gmail.com');
setpref('Internet', 'SMTP_Username', 'manuel.ade@gmail.com');
setpref('Internet', 'SMTP_Password', 'xxxxxx');
setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port', '465');

if x == 'yes'; % Send message
sendmail('manuel.ade@gmail.com', 'Mail from XXX', ...
    ['Hello Mr. Diaz! This is a mesage from XXX!' 10 ...
     'Your simulation has completed succesfully']);
end
 