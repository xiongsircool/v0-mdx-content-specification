import fs from 'fs/promises'
import path from 'path'

// 缓存文件路径
const CACHE_DIR = path.join(process.cwd(), '.cache', 'literature')
const CACHE_FILE = path.join(CACHE_DIR, 'latest-papers.json')
const CACHE_DURATION = 24 * 60 * 60 * 1000 // 24小时缓存

// 外部论文数据接口
export interface ExternalPaper {
  thesisId: string
  paperId: string
  title: string
  titleZh?: string
  summary: string
  summaryZh?: string
  journal?: string
  doi: string
  pmid?: string
  publishDate: string
  pubDate?: string
  authors?: string[]
  keywords?: string[]
  impactFactor?: number
  citations?: number
  volume?: string
  issue?: string
}

interface CachedData {
  papers: ExternalPaper[]
  timestamp: number
  total: number
}

// 确保缓存目录存在
async function ensureCacheDir() {
  try {
    await fs.access(CACHE_DIR)
  } catch {
    await fs.mkdir(CACHE_DIR, { recursive: true })
  }
}

// 从外部API获取论文数据
async function fetchExternalPapers(
  keyword: string = "bioinformatics",
  page: number = 1,
  limit: number = 10
): Promise<{ papers: ExternalPaper[], total: number }> {
  const apiUrl = new URL("https://www.x-mol.net/api/u/paper/search")
  apiUrl.searchParams.append("option", keyword)
  apiUrl.searchParams.append("pageIndex", page.toString())
  apiUrl.searchParams.append("searchSort", "publishDate")

  const response = await fetch(apiUrl.toString(), {
    headers: {
      "Accept": "application/json",
      "User-Agent": "Mozilla/5.0 (compatible; SBC-Bioinformatics/1.0)"
    }
  })

  if (!response.ok) {
    throw new Error(`API request failed: ${response.status} ${response.statusText}`)
  }

  const data = await response.json()
  const papers: ExternalPaper[] = []

  if (data.value?.paperSimpleSearchResult?.pageResults?.results) {
    data.value.paperSimpleSearchResult.pageResults.results.forEach((item: any) => {
      papers.push({
        thesisId: item.thesisId,
        paperId: item.paperId,
        title: item.title,
        titleZh: item.titleZh,
        summary: item.summary,
        summaryZh: item.summaryZh,
        journal: item.journal,
        doi: item.doi,
        pmid: item.pmid,
        publishDate: item.publishDate,
        pubDate: item.pubDate,
        keywords: item.keywordList,
        volume: item.volume,
        issue: item.issue
      })
    })
  }

  return {
    papers,
    total: data.value?.paperSimpleSearchResult?.pageResults?.totalRecord || 0
  }
}

// 检查缓存是否有效
async function isCacheValid(): Promise<boolean> {
  try {
    await ensureCacheDir()
    const content = await fs.readFile(CACHE_FILE, 'utf-8')
    const cached: CachedData = JSON.parse(content)
    return Date.now() - cached.timestamp < CACHE_DURATION
  } catch {
    return false
  }
}

// 获取缓存的论文数据
export async function getCachedPapers(
  forceRefresh: boolean = false
): Promise<{ papers: ExternalPaper[], total: number, fromCache: boolean }> {
  await ensureCacheDir()

  // 检查缓存是否有效
  if (!forceRefresh && await isCacheValid()) {
    try {
      const content = await fs.readFile(CACHE_FILE, 'utf-8')
      const cached: CachedData = JSON.parse(content)
      return {
        papers: cached.papers,
        total: cached.total,
        fromCache: true
      }
    } catch (error) {
      console.error('Error reading cache:', error)
    }
  }

  // 缓存无效或强制刷新，从API获取新数据
  try {
    const { papers, total } = await fetchExternalPapers("bioinformatics", 1, 10)

    // 保存到缓存
    const cachedData: CachedData = {
      papers: papers.slice(0, 6), // 只缓存前6篇用于展示
      total,
      timestamp: Date.now()
    }

    await fs.writeFile(CACHE_FILE, JSON.stringify(cachedData, null, 2))

    return {
      papers: cachedData.papers,
      total: cachedData.total,
      fromCache: false
    }
  } catch (error) {
    console.error('Error fetching external papers:', error)

    // 如果API调用失败，尝试返回旧的缓存
    try {
      const content = await fs.readFile(CACHE_FILE, 'utf-8')
      const cached: CachedData = JSON.parse(content)
      return {
        papers: cached.papers,
        total: cached.total,
        fromCache: true
      }
    } catch {
      // 返回空数组
      return {
        papers: [],
        total: 0,
        fromCache: false
      }
    }
  }
}

// 清除缓存
export async function clearLiteratureCache(): Promise<void> {
  try {
    await fs.unlink(CACHE_FILE)
  } catch (error) {
    // 文件不存在或其他错误，忽略
  }
}

// 手动触发数据更新
export async function updateLiteratureData(): Promise<boolean> {
  try {
    await getCachedPapers(true)
    return true
  } catch (error) {
    console.error('Error updating literature data:', error)
    return false
  }
}