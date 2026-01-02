FROM node:20-alpine

WORKDIR /app

# Install dependencies
COPY package*.json ./
RUN npm install

# Copy project files
COPY . .

# Expose port for dev server
EXPOSE 4321

# Start dev server
CMD ["npm", "run", "dev", "--", "--host"]
