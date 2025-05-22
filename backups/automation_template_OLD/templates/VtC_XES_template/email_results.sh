#!/bin/bash

# Prompt user for email body
read -p "Enter email body (press Enter for default): " body
body=${body:-"Hello,

This is email was auto generated.

Sincerely,
Charles Cardot"}

# Prompt user for subject
read -p "Enter email subject: " subject

# Prompt user for file(s) to attach
read -p "Enter path to file(s) to attach (separated by space, or use 'ALL' for any .png or .pdf in this folder): " files


if [[ $files == "ALL" ]]; then
	out=$(find . -maxdepth 1 \( -name "*.pdf" -o -name "*.png" \))
	files=""
	for file in $out; do
		files+="$file -a "
	done
	files="${files::-3}"
fi

# Prompt user for email recipient(s)
read -p "Enter email recipient(s) (separated by commas): " recipients

# Construct the mail command
mail_command="echo -e \"$body\" | mail -s \"$subject\" -S smtp=smtp://smtp.elasticemail.com:2525 -S smtp-use-starttls -S smtp-auth=plain -S smtp-auth-user=LinuxEmails@linux.com -S smtp-auth-password=230C552E0A26003305B6F0737721A920C6AC -r charleslinux@proton.me -a $files \"$recipients\""

# Print the constructed mail command
echo "The following command will be executed:"
echo "$mail_command"

# Send email
eval "$mail_command"
