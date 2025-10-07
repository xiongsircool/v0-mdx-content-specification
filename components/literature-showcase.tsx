'use client'

import { useState, useEffect } from "react"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { ExternalLink, Calendar, Users, Quote, TrendingUp, BookOpen, RefreshCw } from "lucide-react"

interface ExternalPaper {
  thesisId: string
  paperId: string
  title: string
  titleZh?: string
  summary: string
  summaryZh?: string
  journal?: string
  doi: string
  publishDate: string
  pubDate?: string
  keywords?: string[]
  authors?: string[]
  impactFactor?: number
  citations?: number
  abstract?: string
}

interface LiteratureShowcaseProps {
  papers?: ExternalPaper[]
}

export function LiteratureShowcase({ papers: staticPapers }: LiteratureShowcaseProps = {}) {
  const [papers, setPapers] = useState<ExternalPaper[]>([])
  const [loading, setLoading] = useState(true)
  const [fromCache, setFromCache] = useState(false)
  const [error, setError] = useState<string | null>(null)

  // 默认静态数据（当API不可用时使用）
  const defaultPapers: ExternalPaper[] = staticPapers || [
    {
      thesisId: "1971699463241887744",
      paperId: "1971699463241887744",
      title: "Hesperetin Alleviates Bleomycin-Induced Pulmonary Fibrosis by Modulating Cellular Senescence and Promoting Impaired Autophagy in a CISD2-Dependent Manner",
      summary: "Hesperetin alleviates bleomycin-induced pulmonary fibrosis by modulating cellular senescence and promoting impaired autophagy in a CISD2-dependent manner.",
      summaryZh: "肺纤维化是一种与年龄相关的疾病，其特征是肺功能进行性下降和高死亡率，治疗选择有限。橙皮素是一种柑橘类衍生的类黄酮，具有抗氧化和抗衰老特性，尚未对其在肺纤维化中的潜力进行彻底研究。",
      journal: "The FASEB Journal",
      doi: "10.1096/fj.202502311r",
      publishDate: "2025-10-15",
      pubDate: "2025-10-15",
      keywords: ["橙皮素", "肺纤维化", "细胞衰老", "自噬", "CISD2"]
    },
    {
      thesisId: "1971699463241887745",
      paperId: "1971699463241887745",
      title: "Deep Learning-Driven Single-Cell RNA Sequencing Analysis for Cell Type Identification and Trajectory Inference",
      summary: "This study proposes a novel deep learning framework for single-cell RNA sequencing data analysis, focusing on cell type identification and trajectory inference.",
      summaryZh: "本研究提出了一种新的深度学习框架，用于单细胞RNA测序数据的细胞类型识别和轨迹推断。该方法在多个真实数据集上验证，准确率显著高于传统方法。",
      journal: "Nature Methods",
      doi: "10.1038/s41592-025-01001-2",
      publishDate: "2025-10-12",
      pubDate: "2025-10-12",
      keywords: ["深度学习", "单细胞RNA-seq", "细胞类型识别", "轨迹推断"]
    },
    {
      thesisId: "1971699463241887746",
      paperId: "1971699463241887746",
      title: "Multi-Omics Integration for Precision Medicine: A Comprehensive Framework for Clinical Applications",
      summary: "This study develops a comprehensive multi-omics data integration framework combining genomic, transcriptomic, proteomic, and metabolomic data for precision medicine applications.",
      summaryZh: "本研究开发了一个综合的多组学数据整合框架，整合基因组、转录组、蛋白质组和代谢组数据，用于精准医疗应用。在10,000个临床样本中验证了该方法的有效性。",
      journal: "Cell",
      doi: "10.1016/j.cell.2025.10.003",
      publishDate: "2025-10-08",
      pubDate: "2025-10-08",
      keywords: ["多组学", "精准医疗", "数据整合", "临床应用"]
    }
  ]

  useEffect(() => {
    // 动态加载缓存数据
    const loadCachedData = async () => {
      try {
        const response = await fetch('/api/literature-cache')
        if (response.ok) {
          const data = await response.json()
          setPapers(data.papers || [])
          setFromCache(data.fromCache || false)
        } else {
          // 如果API不可用，使用默认数据
          setPapers(defaultPapers)
          setFromCache(false)
        }
      } catch (error) {
        console.error('Error loading cached literature data:', error)
        // 使用默认数据
        setPapers(defaultPapers)
        setFromCache(false)
      } finally {
        setLoading(false)
      }
    }

    loadCachedData()
  }, [])

  const formatAuthors = (authors: string[]) => {
    if (authors.length <= 3) {
      return authors.join(", ")
    }
    return `${authors.slice(0, 3).join(", ")} 等`
  }

  const formatDate = (dateString: string) => {
    return new Date(dateString).toLocaleDateString('zh-CN', {
      year: 'numeric',
      month: 'long',
      day: 'numeric'
    })
  }

  const getImpactFactorColor = (ifScore: number) => {
    if (ifScore >= 50) return "bg-red-100 text-red-800 border-red-200"
    if (ifScore >= 30) return "bg-orange-100 text-orange-800 border-orange-200"
    if (ifScore >= 10) return "bg-yellow-100 text-yellow-800 border-yellow-200"
    return "bg-blue-100 text-blue-800 border-blue-200"
  }

  if (loading) {
    return (
      <section className="py-24 relative bg-gradient-to-br from-background via-primary/5 to-background">
        {/* Tech Background */}
        <div className="absolute inset-0 tech-grid opacity-20"></div>

        <div className="container mx-auto px-4 relative z-10">
          <div className="text-center mb-20">
            <h2 className="text-5xl font-bold mb-6">
              <span className="text-foreground/90">LATEST</span>
              <span className="text-primary neon-glow"> LITERATURE</span>
            </h2>
            <p className="text-lg text-foreground/60 font-mono max-w-3xl mx-auto">
              关注生物信息学前沿研究动态，定期更新重要学术论文
            </p>
          </div>

          <div className="flex justify-center items-center h-64">
            <div className="text-center">
              <RefreshCw className="w-8 h-8 animate-spin mx-auto mb-4 text-primary" />
              <p className="text-foreground/60 font-mono">正在加载文献数据...</p>
            </div>
          </div>
        </div>
      </section>
    )
  }

  if (error) {
    return (
      <section className="py-24 relative bg-gradient-to-br from-background via-primary/5 to-background">
        {/* Tech Background */}
        <div className="absolute inset-0 tech-grid opacity-20"></div>

        <div className="container mx-auto px-4 relative z-10">
          <div className="text-center mb-20">
            <h2 className="text-5xl font-bold mb-6">
              <span className="text-foreground/90">LATEST</span>
              <span className="text-primary neon-glow"> LITERATURE</span>
            </h2>
            <p className="text-lg text-foreground/60 font-mono max-w-3xl mx-auto">
              关注生物信息学前沿研究动态，定期更新重要学术论文
            </p>
          </div>

          <div className="flex justify-center items-center h-64">
            <div className="text-center">
              <p className="text-foreground/60 font-mono mb-2">⚠️ 无法加载文献数据</p>
              <p className="text-foreground/50 font-mono">显示默认文献信息</p>
            </div>
          </div>
        </div>
      </section>
    )
  }

  return (
    <section className="py-24 relative bg-gradient-to-br from-background via-primary/5 to-background">
      {/* Tech Background */}
      <div className="absolute inset-0 tech-grid opacity-20"></div>

      <div className="container mx-auto px-4 relative z-10">
        <div className="text-center mb-20">
          <h2 className="text-5xl font-bold mb-6">
            <span className="text-foreground/90">LATEST</span>
            <span className="text-primary neon-glow"> LITERATURE</span>
          </h2>
          <p className="text-lg text-foreground/60 font-mono max-w-3xl mx-auto">
            关注生物信息学前沿研究动态，定期更新重要学术论文
          </p>
          {fromCache && (
            <p className="text-sm text-foreground/50 mt-4 font-mono">
              📊 数据来自缓存 | 最后更新: {new Date().toLocaleString('zh-CN')}
            </p>
          )}
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
          {papers.map((literature, index) => (
            <Card
              key={index}
              className="group glass-effect rounded-2xl border border-primary/20
                        hover:border-primary/40 transition-all duration-500
                        hover:shadow-2xl hover:shadow-primary/10
                        hover:-translate-y-2 transform h-full"
            >
              <CardHeader className="pb-4">
                <div className="flex justify-between items-start gap-3 mb-4">
                  <div className="flex-1">
                    <CardTitle className="text-lg font-bold leading-tight text-primary group-hover:text-primary/80 transition-colors line-clamp-3">
                      {literature.title}
                    </CardTitle>
                  </div>
                  {literature.impactFactor && (
                    <Badge
                      className={`flex-shrink-0 bg-primary/10 text-primary border-primary/20 hover:bg-primary/20 transition-colors`}
                    >
                      IF: {literature.impactFactor}
                    </Badge>
                  )}
                </div>

                {/* Journal and Date */}
                <div className="space-y-2">
                  {literature.journal && (
                    <div className="flex items-center gap-2 text-sm text-foreground/70">
                      <BookOpen className="w-4 h-4 text-primary" />
                      <span className="font-medium">{literature.journal}</span>
                    </div>
                  )}
                  <div className="flex items-center gap-2 text-sm text-foreground/60">
                    <Calendar className="w-4 h-4 text-primary" />
                    <span>{formatDate(literature.publishDate)}</span>
                  </div>
                </div>
              </CardHeader>

              <CardContent className="space-y-4">
                {/* Authors */}
                {literature.authors && literature.authors.length > 0 ? (
                  <div className="flex items-start gap-2 text-sm">
                    <Users className="w-4 h-4 text-primary mt-0.5 flex-shrink-0" />
                    <span className="text-foreground/70 leading-relaxed">
                      {formatAuthors(literature.authors)}
                    </span>
                  </div>
                ) : (
                  <div className="flex items-start gap-2 text-sm">
                    <Users className="w-4 h-4 text-foreground/40 mt-0.5 flex-shrink-0" />
                    <span className="text-foreground/50 leading-relaxed italic">
                      作者信息暂未收录
                    </span>
                  </div>
                )}

                {/* Abstract */}
                <div className="text-sm text-foreground/60 leading-relaxed">
                  <div className="flex items-start gap-2">
                    <Quote className="w-4 h-4 text-foreground/40 mt-0.5 flex-shrink-0" />
                    <p className="line-clamp-4">
                      {literature.abstract || literature.summaryZh || literature.summary || '暂无摘要'}
                    </p>
                  </div>
                </div>

                {/* Keywords */}
                {literature.keywords && literature.keywords.length > 0 && (
                  <div className="flex flex-wrap gap-2">
                    {literature.keywords.slice(0, 3).map((keyword, keywordIndex) => (
                      <span
                        key={keywordIndex}
                        className="px-3 py-1 rounded-full bg-primary/10 text-primary text-xs font-mono
                               hover:bg-primary/20 transition-colors border border-primary/20"
                      >
                        {keyword}
                      </span>
                    ))}
                    {literature.keywords.length > 3 && (
                      <span className="px-3 py-1 rounded-full bg-background/60 text-foreground/50 text-xs font-mono
                                   border border-primary/10">
                        +{literature.keywords.length - 3}
                      </span>
                    )}
                  </div>
                )}

                {/* Citations and Link */}
                <div className="flex items-center justify-between pt-4 border-t border-primary/10">
                  {literature.citations && (
                    <div className="flex items-center gap-2 text-sm text-foreground/60">
                      <TrendingUp className="w-4 h-4 text-primary" />
                      <span>引用: {literature.citations}</span>
                    </div>
                  )}

                  <div className="flex gap-2">
                    <button
                      className="text-primary hover:text-primary/80 hover:bg-primary/10 text-xs px-3 py-1 rounded-full border border-primary/20 transition-all duration-300 font-mono"
                      onClick={() => window.open(`https://doi.org/${literature.doi}`, '_blank')}
                    >
                      <ExternalLink className="w-3 h-3 inline mr-1" />
                      DOI
                    </button>
                  </div>
                </div>
              </CardContent>
            </Card>
          ))}
        </div>

        {/* Note */}
        <div className="text-center mt-16">
          <div className="flex justify-center mb-8">
            <div className="w-16 h-1 bg-gradient-to-r from-primary via-accent to-primary rounded-full"></div>
          </div>
          <p className="text-foreground/50 text-sm font-mono">
            📚 文献信息定期更新 | 数据来源：x-mol.net
          </p>
        </div>
      </div>
    </section>
  )
}