#!/bin/bash

# Prompt user for subject
read -p "Enter email subject: " subject

# Prompt user for email body
read -p "Enter email body (press Enter for default): " body
body=${body:-"<p>Hello,<br><br>This email was auto-generated.<br><br>Sincerely,<br>Charles Cardot</p>"}

# Prompt user for file(s) to attach
read -p "Enter path to file(s) to attach (separated by space, or use 'ALL' for any .png or .pdf in this folder): " files

if [[ $files == "ALL" ]]; then
  out=$(find . -maxdepth 1 \( -name "*.pdf" -o -name "*.png" \))
  files=""
  for file in $out; do
      files+="$file "
  done
fi

# Initialize attachments array
attachments=""

# Process each file
for file in $files; do
  # Encode the file to Base64 and store it in a temporary location
  base64 -w 0 "$file" > /tmp/file_base64.txt
  BASE64_CONTENT=$(cat /tmp/file_base64.txt)
  filename=$(basename "$file")

  # Append attachment information to the attachments array
  attachments+="{ \"filename\": \"$filename\", \"content\": \"$BASE64_CONTENT\" },"
done

# Remove trailing comma from attachments array
attachments="${attachments%,}"

# Prompt user for multiple comma-separated email recipients
read -p "Enter email recipients (comma-separated), (WARNING: Mailsender does not support more than one recipient): " recipients

# Format recipients into JSON array
recipients_json="["
IFS=',' read -r -a recipient_array <<< "$recipients"
for recipient in "${recipient_array[@]}"; do
  recipient=$(echo "$recipient" | xargs)  # Trim leading and trailing whitespace
  if [[ -n "$recipient" ]]; then
    recipients_json+="{ \"email\": \"$recipient\" },"
  fi
done
recipients_json="${recipients_json%,}]"

# Create the JSON request body and store it in a temporary file
echo "{
  \"from\": { \"email\": \"MS_RddIyf@trial-k68zxl258qelj905.mlsender.net\" },
  \"to\": $recipients_json,
  \"subject\": \"$subject\",
  \"html\": \"$body\",
  \"attachments\": [$attachments]
}" > /tmp/request_body.json

# Construct the curl command to send the email via MailerSend
mail_command="curl -X POST https://api.mailersend.com/v1/email \
  -H 'Content-Type: application/json' \
  -H 'X-Requested-With: XMLHttpRequest' \
  -H 'Authorization: Bearer mlsn.e3ac2bdf065bf2b0f805b36a80e28df0165a394baf851061ce49378412feda4d' \
  --data-binary @/tmp/request_body.json"

# Print the constructed mail command
echo "The following command will be executed:"
echo "$mail_command"

# Send email
eval "$mail_command"

# Clean up temporary files
rm /tmp/file_base64.txt
rm /tmp/request_body.json
