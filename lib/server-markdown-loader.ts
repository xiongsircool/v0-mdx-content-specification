import matter from 'gray-matter'
import { marked } from 'marked'
import hljs from 'highlight.js'

// 只在服务端导入fs
let fs: any, path: any
if (typeof window === 'undefined') {
  fs = require('fs')
  path = require('path')
}

// 配置 marked 以支持代码高亮
marked.setOptions({
  highlight: function(code, lang) {
    if (lang && hljs.getLanguage(lang)) {
      try {
        return hljs.highlight(code, { language: lang }).value
      } catch (err) {
        console.error('Highlight.js error:', err)
      }
    }
    return hljs.highlightAuto(code).value
  },
  breaks: true,
  gfm: true
})

export interface MarkdownMetadata {
  title: string
  publishedAt?: string
  updatedAt?: string
  date?: string
  time?: string
  location?: string
  eventType?: string
  status?: string
  summary?: string
  excerpt?: string
  tags?: string[]
  authors?: string[]
  coverImage?: string
  type?: string
  level?: string
  track?: string
  duration?: string
  prerequisites?: string[]
  externalUrl?: string
  priority?: number
  expiresAt?: string
  attachments?: Array<{name: string, url: string}>
  speaker?: string
  capacity?: number
  registered?: number
  meetingUrl?: string
  materials?: string[]
  [key: string]: any
}

export interface MarkdownFile {
  content: string
  metadata: MarkdownMetadata
  htmlContent: string
  slug: string
  url: string
  filePath: string
}

export interface ContentItem {
  title: string
  url: string
  slug: string
  summary?: string
  excerpt?: string
  publishedAt?: string
  updatedAt?: string
  date?: string
  tags?: string[]
  authors?: string[]
  type?: string
  level?: string
  track?: string
  coverImage?: string
  priority?: number
  expiresAt?: string
  isExpired?: boolean
  [key: string]: any
}

// 解析单个markdown文件
export function parseMarkdownFile(filePath: string, baseUrl: string): MarkdownFile {
  const fileContent = fs.readFileSync(filePath, 'utf8')
  const { data, content } = matter(fileContent)

  // 生成slug和url
  const relativePath = filePath.replace(process.cwd() + '/content/', '')
  const slug = relativePath.replace(/\.(md|markdown)$/, '')
  const url = `${baseUrl}/${slug}`

  // 转换markdown为HTML
  const htmlContent = marked(content)

  return {
    content,
    metadata: data as MarkdownMetadata,
    htmlContent,
    slug,
    url,
    filePath
  }
}

// 加载指定目录下的所有markdown文件
export function loadMarkdownFiles(directory: string, baseUrl: string): MarkdownFile[] {
  const fullPath = path.join(process.cwd(), 'content', directory)

  if (!fs.existsSync(fullPath)) {
    return []
  }

  const files = fs.readdirSync(fullPath)
  const markdownFiles: MarkdownFile[] = []

  files.forEach(file => {
    const filePath = path.join(fullPath, file)
    const stat = fs.statSync(filePath)

    if (stat.isFile() && (file.endsWith('.md') || file.endsWith('.markdown'))) {
      try {
        const parsedFile = parseMarkdownFile(filePath, baseUrl)
        markdownFiles.push(parsedFile)
      } catch (error) {
        console.error(`Error parsing ${filePath}:`, error)
      }
    }
  })

  return markdownFiles
}

// 获取公告列表
export function getAnnouncements(): ContentItem[] {
  const files = loadMarkdownFiles('announcements', '')

  return files.map(file => ({
    title: file.metadata.title,
    url: file.url,
    slug: file.slug,
    summary: file.metadata.summary,
    publishedAt: file.metadata.publishedAt,
    tags: file.metadata.tags,
    priority: file.metadata.priority || 0,
    expiresAt: file.metadata.expiresAt,
    isExpired: file.metadata.expiresAt ? new Date(file.metadata.expiresAt) < new Date() : false
  })).sort((a, b) => {
    // 按优先级和发布日期排序
    if (a.priority !== b.priority) {
      return (b.priority || 0) - (a.priority || 0)
    }
    return new Date(b.publishedAt || '').getTime() - new Date(a.publishedAt || '').getTime()
  })
}

// 获取文章列表
export function getPosts(): ContentItem[] {
  const files = loadMarkdownFiles('posts', '')

  return files.map(file => ({
    title: file.metadata.title,
    url: file.url,
    slug: file.slug,
    excerpt: file.metadata.excerpt,
    publishedAt: file.metadata.publishedAt,
    updatedAt: file.metadata.updatedAt,
    tags: file.metadata.tags,
    authors: file.metadata.authors,
    coverImage: file.metadata.coverImage
  })).sort((a, b) => {
    return new Date(b.publishedAt || '').getTime() - new Date(a.publishedAt || '').getTime()
  })
}

// 获取学习资源列表
export function getLearnResources(): ContentItem[] {
  const files = loadMarkdownFiles('learn/resources', '')

  return files.map(file => ({
    title: file.metadata.title,
    url: file.url,
    slug: file.slug,
    type: file.metadata.type,
    level: file.metadata.level,
    track: file.metadata.track,
    duration: file.metadata.duration,
    prerequisites: file.metadata.prerequisites,
    externalUrl: file.metadata.externalUrl,
    updatedAt: file.metadata.updatedAt
  })).sort((a, b) => {
    return new Date(b.updatedAt || '').getTime() - new Date(b.updatedAt || '').getTime()
  })
}

// 获取会议列表
export function getMeetings(): ContentItem[] {
  const files = loadMarkdownFiles('meetings', '')

  return files.map(file => ({
    title: file.metadata.title,
    url: file.url,
    slug: file.slug,
    date: file.metadata.date,
    time: file.metadata.time,
    location: file.metadata.location,
    eventType: file.metadata.eventType,
    status: file.metadata.status,
    summary: file.metadata.summary,
    tags: file.metadata.tags,
    speaker: file.metadata.speaker,
    capacity: file.metadata.capacity,
    registered: file.metadata.registered,
    meetingUrl: file.metadata.meetingUrl,
    materials: file.metadata.materials
  })).sort((a, b) => {
    return new Date(a.date || '').getTime() - new Date(b.date || '').getTime()
  })
}

// 根据slug获取单个文件
export function getFileBySlug(type: string, slug: string): MarkdownFile | null {
  const possiblePaths = [
    `content/${type}/${slug}.md`,
    `content/${type}/${slug}.markdown`
  ]

  for (const filePath of possiblePaths) {
    const fullPath = path.join(process.cwd(), filePath)
    if (fs.existsSync(fullPath)) {
      try {
        return parseMarkdownFile(fullPath, '')
      } catch (error) {
        console.error(`Error parsing ${fullPath}:`, error)
      }
    }
  }

  return null
}

// 搜索功能
export function searchContent(query: string): ContentItem[] {
  const allContent = [
    ...getAnnouncements(),
    ...getPosts(),
    ...getLearnResources(),
    ...getMeetings()
  ]

  const lowercaseQuery = query.toLowerCase()

  return allContent.filter(item => {
    return (
      item.title.toLowerCase().includes(lowercaseQuery) ||
      item.summary?.toLowerCase().includes(lowercaseQuery) ||
      item.excerpt?.toLowerCase().includes(lowercaseQuery) ||
      item.tags?.some(tag => tag.toLowerCase().includes(lowercaseQuery))
    )
  })
}